import esutil, fitsio
import healpy as hp
import numpy as np
from catalog import Catalog,Entry
from utilities import TOTAL_SQDEG, SEC_PER_DEG, astro_to_sphere, calc_theta_i, apply_errormodels
from scipy.special import erf
#from cluster import Cluster

class Mask(object):
    """
    A super-class for (pixelized) footprint masks

    This should not be instantiated directly (yet).

    parameters
    ----------
    confstr: Config object
       configuration
    """

    # note: not sure how to organize this.
    #   We need a routine that looks at the mask_mode and instantiates
    #   the correct type.  How is this typically done?

    def __init__(self, confstr):
        try:
            self.read_maskgals(confstr.maskgalfile)
        except:
            # this could throw a ValueError or AttributeError
            self.gen_maskgals()

    def calc_radmask(self, *args, **kwargs): pass
    def read_maskgals(self, maskgalfile):
        self.maskgals = Catalog.from_fits_file(maskgalfile)
    def gen_maskgals(self):
        # this needs to be written to generate maskgals if not from file
        # Tom-where would we generate them from?
        pass

    def set_radmask(self, cluster, mpcscale):
        """
        Assign mask (0/1) values to maskgals for a given cluster

        parameters
        ----------
        cluster: Cluster object
        mpcscale: float
            scaling to go from mpc to degrees (check units) at cluster redshift

        results
        -------
        sets maskgals.mark

        """

        # note this probably can be in the superclass, no?
        ras = cluster.ra + self.maskgals.x/(mpcscale*SEC_PER_DEG)/np.cos(np.radians(cluster.dec))
        decs = cluster.dec + self.maskgals.y/(mpcscale*SEC_PER_DEG)
        self.maskgals.mark = self.compute_radmask(ras,decs)

    def calc_maskcorr(self, mstar, maxmag, limmag):
        """
        Obtain mask correction c parameters. From calclambda_chisq_calc_maskcorr.pro

        parameters
        ----------
        maskgals : Object holding mask galaxy parameters
        mstar    :
        maxmag   : Maximum magnitude
        limmag   : Limiting Magnitude
        confstr  : Configuration object
                    containing configuration info

        returns
        -------
        cpars

        """

        mag_in = self.maskgals.m + mstar
        self.maskgals.refmag = mag_in

        if self.maskgals.limmag[0] > 0.0:
            mag, mag_err = apply_errormodels(self.maskgals, mag_in)

            self.maskgals.refmag_obs = mag
            self.maskgals.refmag_obs_err = mag_err
        else:
            mag = mag_in
            mag_err = 0*mag_in
            raise ValueError('Survey limiting magnitude <= 0!')
            #Raise error here as this would lead to divide by zero if called.

        #extract object for testing
        #fitsio.write('test_data.fits', self.maskgals._ndarray)

        if (self.maskgals.w[0] < 0) or (self.maskgals.w[0] == 0 and 
            np.amax(self.maskgals.m50) == 0):
            theta_i = calc_theta_i(mag, mag_err, maxmag, limmag)
        elif (self.maskgals.w[0] == 0):
            theta_i = calc_theta_i(mag, mag_err, maxmag, self.maskgals.m50)
        else:
            raise Exception('Unsupported mode!')

        p_det = theta_i*self.maskgals.mark
        np.set_printoptions(threshold=np.nan)
        #print self.maskgals.mark
        c = 1 - np.dot(p_det, self.maskgals.theta_r) / self.maskgals.nin[0]

        cpars = np.polyfit(self.maskgals.radbins[0], c, 3)

        return cpars

class HPMask(Mask):
    """
    A class to use a healpix mask (mask_mode == 3)

    parameters
    ----------
    confstr: Config object
        Configuration object with maskfile

    """

    def __init__(self, confstr):
        # record for posterity
        self.maskfile = confstr.maskfile
        maskinfo, hdr = fitsio.read(confstr.maskfile, ext=1, header=True)
        # maskinfo converted to a catalog (array of Entrys)
        maskinfo = Catalog(maskinfo)
        nlim, nside, nest = maskinfo.hpix.size, hdr['NSIDE'], hdr['NEST']
        hpix_ring = maskinfo.hpix if nest != 1 else hp.nest2ring(nside, maskinfo.hpix)
        muse = np.arange(nlim)

        # if we have a sub-region of the sky, cut down the mask to save memory
        if confstr.hpix > 0:
            border = confstr.border + hp.nside2resol(nside)
            theta, phi = hp.pix2ang(confstr.nside, confstr.hpix)
            radius = np.sqrt(2) * (hp.nside2resol(confstr.nside)/2. + border)
            pixint = hp.query_disc(nside, hp.ang2vec(theta, phi), 
                                        np.radians(radius), inclusive=False)
            muse, = esutil.numpy_util.match(hpix_ring, pixint)

        offset, ntot = np.min(hpix_ring)-1, np.max(hpix_ring)-np.min(hpix_ring)+3
        self.nside = nside
        self.offset = offset
        self.npix = ntot

        #ntot = np.max(hpix_ring) - np.min(hpix_ring) + 3
        self.fracgood = np.zeros(ntot,dtype='f4')

        # check if we have a fracgood in the input maskinfo
        try:
            self.fracgood_float = 1
            self.fracgood[hpix_ring-offset] = maskinfo[muse].fracgood
        except AttributeError:
            self.fracgood_float = 0
            self.fracgood[hpix_ring-offset] = 1
        super(HPMask, self).__init__(confstr)

    def compute_radmask(self, ra, dec):
        """
        Determine if a given set of ra/dec points are in or out of mask

        parameters
        ----------
        ra: array of doubles
        dec: array of doubles

        returns
        -------
        radmask: array of booleans

        """
        _ra  = np.atleast_1d(ra)
        _dec = np.atleast_1d(dec)

        if (_ra.size != _dec.size):
            raise ValueError("ra, dec must be same length")

        theta, phi = astro_to_sphere(_ra, _dec)
        ipring = hp.ang2pix(self.nside, theta, phi)
        ipring_offset = np.clip(ipring - self.offset, 0, self.npix-1)
        ref = 0 if self.fracgood_float == 0 else np.random.rand(_ra.size)
        radmask = np.zeros(_ra.size, dtype=np.bool_)
        radmask[np.where(self.fracgood[ipring_offset] > ref)] = True
        return radmask
