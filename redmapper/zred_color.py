from __future__ import division, absolute_import, print_function
from past.builtins import xrange

import numpy as np
import esutil
import scipy.integrate
import scipy.interpolate
import copy

from .galaxy import GalaxyCatalog
from .utilities import interpol

class ZredColor(object):
    """
    """

    def __init__(self, zredstr, sigint=0.001, do_correction=True, adaptive=True,
                 use_photoerr=True, zrange=None):
        self.zredstr = zredstr

        self.sigint = sigint
        self.do_correction = do_correction
        self.use_photoerr = use_photoerr
        self.zrange = zrange

        self.nz = self.zredstr.z.size - 1
        self.notextrap, = np.where(~self.zredstr.extrapolated)

        if self.zrange is None:
            self.zbinstart = 0
            self.zbinstop = self.nz - 1
        else:
            u, = np.where((self.zredstr.z > self.zrange[0]) &
                          (self.zredstr.z < self.zrange[1]))
            self.zbinstart = u[0]
            self.zbinstop = u[-1]

    def compute_zreds(self, galaxies):
        """
        """

        for galaxy in galaxies:
            self.compute_zred(galaxy, no_corrections=True)

        if self.do_correction:
            # Bulk processing
            olddzs = np.zeros(galaxies.size)
            dzs = np.zeros_like(olddzs)
            iteration = 0

            pivotmags = interpol(self.zredstr.pivotmag, self.zredstr.z, galaxies.zred_uncorr)

            while (iteration < 5):
                olddzs[:] = dzs
                dzs[:] = interpol(self.zredstr.corr, self.zredstr.z, galaxies.zred_uncorr + olddzs) + (galaxies.refmag - pivotmags) * interpol(self.zredstr.corr_slope, self.zredstr.z, galaxies.zred_uncorr + olddzs)
                iteration += 1

            galaxies.zred = galaxies.zred_uncorr + dzs
            galaxies.zred_e = galaxies.zred_uncorr_e * interpol(self.zredstr.corr_r, self.zredstr.z, galaxies.zred)

            dz2s = interpol(self.zredstr.corr2, self.zredstr.z, galaxies.zred_uncorr)
            r2s = interpol(self.zredstr.corr2_r, self.zredstr.z, galaxies.zred_uncorr)

            galaxies.zred2 = galaxies.zred_uncorr + dz2s
            galaxies.zred2_e = galaxies.zred_uncorr_e * r2s


    def compute_zred(self, galaxy, no_corrections=False):
        """
        """

        lndist = np.zeros(self.nz) - 1e12
        chisq = np.zeros(self.nz) + 1e12

        zbins_limited, = np.where((galaxy.refmag < self.zredstr.maxrefmag) &
                                  (galaxy.refmag > self.zredstr.minrefmag) &
                                  (self.zredstr.z < 100.0))

        if zbins_limited.size < 2:
            self._reset_bad_values(galaxy)
            return

        neighbors = 10
        zbins = np.arange(*np.clip([zbins_limited[0] - neighbors,
                                    zbins_limited[-1] + neighbors],
                                   0, self.nz))

        lndist[zbins], chisq[zbins] = self._calculate_lndist(galaxy, zbins)

        # move from log space to regular space
        maxlndist = np.max(lndist[zbins_limited])
        dist = np.zeros_like(lndist)
        with np.errstate(invalid='ignore', over='ignore'):
            dist[zbins_limited] = np.exp(lndist[zbins_limited] - maxlndist)

        # fix infinities and NaNs
        bad, = np.where(~np.isfinite(dist))
        dist[bad] = 0.0

        # take the maximum where not extrapolated
        ind_temp = np.argmax(dist[self.notextrap])
        ind = self.notextrap[ind_temp]

        calcinds, = np.where(dist > 1e-5)

        if calcinds.size >= 3:
            tdist = scipy.integrate.trapz(dist[calcinds], self.zredstr.z[calcinds])
            zred_temp = scipy.integrate.trapz(dist[calcinds] * self.zredstr.z[calcinds],
                                              self.zredstr.z[calcinds]) / tdist
            zred_e = scipy.integrate.trapz(dist[calcinds] * self.zredstr.z[calcinds]**2.,
                                           self.zredstr.z[calcinds]) / tdist - zred_temp**2.
        else:
            tdist = np.sum(dist[calcinds])
            zred_temp = np.sum(dist[calcinds] * self.zredstr.z[calcinds]) / tdist
            zred_e = np.sum(dist[calcinds] * self.zredstr.z[calcinds]**2.) / tdist - zred_temp**2.

        if zred_e < 0.0:
            zred_e = 1.0
        else:
            zred_e = np.sqrt(zred_e)

        zred_e = zred_e if zred_e > 0.005 else 0.005

        # Now fit a parabola to get the perfect zred
        ind_temp = np.argmax(dist[self.notextrap])
        ind = self.notextrap[ind_temp]

        zred = zred_temp.copy()

        neighbors = 2
        use, = np.where(lndist > -1e10)
        if use.size >= neighbors * 2 + 1:
            minindex, maxindex = np.clip([ind - neighbors, ind + neighbors], use[0], use[-1])

            if ((maxindex - minindex + 1) >= 5):
                X = np.zeros((maxindex - minindex + 1, 3))
                X[:, 1] = self.zredstr.z[minindex:maxindex + 1]
                X[:, 0] = X[:, 1] * X[:, 1]
                X[:, 2] = 1
                y = lndist[minindex: maxindex + 1]

                fit = np.matmul(np.matmul(np.linalg.inv(np.matmul(X.T, X)), X.T), y)

                if fit[0] < 0.0:
                    ztry = -fit[1] / (2.0 * fit[0])
                    # Don't let it move to far, or it's a bad fit
                    if (np.abs(ztry - zred) < 2.0*zred_e):
                        zred = ztry

        # And compute values at the real zred peak
        x = (self.zredstr.z - zred) / zred_e
        newdist = np.exp(-0.5 * x * x)

        bad, = np.where((lndist < -1e10) | (~np.isfinite(lndist)))
        newdist[bad] = 0.0
        lndist[bad] = -1e11

        if calcinds.size >= 3:
            # Note there maybe should be a distcorr here, but this is not
            #  actually computed in the IDL code (bug?)
            lkhd = scipy.integrate.trapz(newdist[calcinds] * (lndist[calcinds]), self.zredstr.z[calcinds]) / scipy.integrate.trapz(newdist[calcinds], self.zredstr.z[calcinds])
        else:
            lkhd = np.sum(newdist[calcinds] * lndist[calcinds]) / np.sum(newdist[calcinds])

        # Get chisq at the closest bin position
        zbin = np.argmin(np.abs(zred - self.zredstr.z))
        chisq = chisq[zbin]

        if not np.isfinite(lkhd):
            self._reset_bad_values(galaxy)
            return

        # And apply the corrections
        zred2 = np.zeros(1) + zred
        zred2_e = np.zeros(1) + zred_e
        zred_uncorr = np.zeros(1) + zred
        zred_uncorr_e = np.zeros(1) + zred_e

        if self.do_correction and not no_corrections:
            olddz = -1.0
            dz = 0.0
            iteration = 0

            pivotmag = interpol(self.zredstr.pivotmag, self.zredstr.z, zred)

            while np.abs(olddz - dz) > 1e-3 and iteration < 10:
                olddz = copy.copy(dz)
                dz = interpol(self.zredstr.corr, self.zredstr.z, zred + olddz) + (galaxy.refmag - pivotmag) * interpol(self.zredstr.corr_slope, self.zredstr.z, zred + olddz)
                iteration += 1

            zred = zred + dz

            # evaluate error correction at "z_true"
            zred_e *= interpol(self.zredstr.corr_r, self.zredstr.z, zred)

            # And the zred2 correction
            dz = interpol(self.zredstr.corr2, self.zredstr.z, zred2) + (galaxy.refmag - pivotmag) * (interpol(self.zredstr.corr2_slope, self.zredstr.z, zred2))
            # this is evaluated at zred0
            r2 = interpol(self.zredstr.corr2_r, self.zredstr.z, zred2)

            zred2 += dz
            zred2_e *= r2

        # Finally store the values

        galaxy.zred = zred
        galaxy.zred_e = zred_e
        galaxy.zred2 = zred2
        galaxy.zred2_e = zred2_e
        galaxy.zred_uncorr = zred_uncorr
        galaxy.zred_uncorr_e = zred_uncorr_e
        galaxy.chisq = chisq
        galaxy.lkhd = lkhd

        # and we're done

    def _calculate_lndist(self, galaxy, zbins):
        """
        """

        # Note we need to deal with photoerr...
        if zbins.size > 1:
            # we have many bins...
            chisq = self.zredstr.calculate_chisq_redshifts(galaxy, zbins, z_is_index=True, calc_lkhd=False)
        else:
            # we have a single bin... hack this
            chisq = self.zredstr.calculate_chisq(galaxy, np.array([zbins[0], zbins[0]]), z_is_index=True, calc_lkhd=False)[0]

        lndist = -0.5 * chisq

        #with np.errstate(invalid='ignore'):
        lndistcorr = np.log((10.**(0.4 * (self.zredstr.alpha + 1.0) *
                                   (self.zredstr._mstar[zbins] - galaxy.refmag)) *
                             np.exp(-10.**(0.4 * (self.zredstr._mstar[zbins] - galaxy.refmag)))) *
                            self.zredstr.volume_factor[zbins])

        lndist += lndistcorr

        bad, = np.where(~np.isfinite(lndist))
        lndist[bad] = -1e11

        return (lndist, chisq)

    def _reset_bad_values(self, galaxy):
        """
        """

        galaxy.lkhd = -1000.0
        galaxy.zred = -1.0
        galaxy.zred_e = -1.0
        galaxy.zred2 = -1.0
        galaxy.zred2_e = -1.0
        galaxy.zred_uncorr = -1.0
        galaxy.zred_uncorr_e = -1.0
        galaxy.chisq = -1.0

