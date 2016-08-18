import unittest
import numpy.testing as testing
import numpy as np
import fitsio

from esutil.cosmology import Cosmo
from redmapper.cluster import Cluster
from redmapper.galaxy import GalaxyCatalog
from redmapper.background import Background
from redmapper.redsequence import RedSequenceColorPar


class BackgroundStub(Background):

    def __init__(self, filename):
        obkg = Entry.from_fits_file(filename)
        self.refmagbins = obkg.refmagbins
        self.chisqbins = obkg.chisqbins
        self.lnchisqbins = obkg.lnchisqbins
        self.zbins = obkg.zbins
        self.sigma_g = obkg.sigma_g
        self.sigma_lng = obkg.sigma_lng


class ClusterFiltersTestCase(unittest.TestCase):

    def test_nfw_filter(self):
        test_indices = np.array([46, 38,  1,  1, 11, 24, 25, 16])
        py_nfw = self._calc_radial_profile()[test_indices]
        idl_nfw = np.array([0.29360449, 0.14824243, 0.14721203, 0.14721203, 
                            0.23459411, 0.31615007, 0.29307860, 0.29737136])
        testing.assert_almost_equal(py_nfw, idl_nfw)

    def test_lum_filter(self):
        test_indices = np.array([47, 19,  0, 30, 22, 48, 34, 19])
        zred_file_name = 'test_dr8_pars.fit'
        zredstr = RedSequenceColorPar(file_path + '/' + zred_file_name)
        py_lum = self._calc_luminosity(zredstr, maxmag)
        idl_lum = np.array([])
        testing.assert_almost_equal(py_lum, idl_lum)

    def test_bkg_filter(self):
        test_indices = np.array([29, 16, 27,  5, 38, 35, 25, 44])
        bkg_file_name = 'test_bkg.fit'
        bkg = BackgroundStub(self.file_path + '/' + bkg_file_name)
        py_bkg = self._calc_bkg_density(bkg, Cosmo())
        idl_bkg = np.array([])
        testing.assert_almost_equal(py_bkg, idl_bkg)

    def setUp(self):
        self.cluster = Cluster(np.empty(0))
        self.file_path, file_name = 'data', 'test_cluster_members.fit'
        self.cluster.members = GalaxyCatalog.from_fits_file(file_path 
                                                            + '/' + file_name)
        self.cluster.z = cluster.members.z[0]

        
class ClusterMembersTestCase(unittest.TestCase):

    def test_member_finding(self): pass
    def test_richness(self): pass


if __name__=='__main__':
    unittest.main()
