from __future__ import division, absolute_import, print_function
from past.builtins import xrange

import matplotlib
matplotlib.use('Agg')

import unittest
import os
import shutil
import numpy.testing as testing
import numpy as np
import fitsio
import tempfile
from numpy import random
import glob

from redmapper import Configuration
from redmapper.redmagic import RedmagicParameterFitter, RedmagicCalibrator
from redmapper.redmagic import RunRedmagicTask
from redmapper import RedSequenceColorPar
from redmapper import GalaxyCatalog
from redmapper import Catalog
from redmapper import VolumeLimitMask

class RedmagicCalTestCase(unittest.TestCase):
    def test_redmagic_fitter(self):
        np.random.seed(12345)

        file_path = 'data_for_tests/redmagic_test'

        # Read in the red-sequence parametrization

        # Read in the input data for comparison (no afterburner)
        calstr = fitsio.read(os.path.join(file_path, 'rcal_str_preab.fit'), ext=1, lower=True)

        # Read in the input data for testing (no afterburner)
        calstr2 = fitsio.read(os.path.join(file_path, 'rcal_str2.fit'), ext=1, lower=True)

        # Make a zred structure for mstar...
        config = Configuration(os.path.join('data_for_tests', 'testconfig.yaml'))
        zredstr = RedSequenceColorPar(None, config=config)

        # Set up the fitter...
        #randomn = np.random.normal(size=calstr2['z'][0, :].size)
        # Old IDL code did not sample for the selection, I think this was wrong
        randomn = np.zeros(calstr2['z'][0, :].size)

        rmfitter = RedmagicParameterFitter(calstr['nodes'][0, :], calstr['corrnodes'][0, :],
                                           calstr2['z'][0, :], calstr2['z_err'][0, :],
                                           calstr2['chisq'][0, :], calstr2['mstar'][0, :],
                                           calstr2['zcal'][0, :], calstr2['zcal_e'][0, :],
                                           calstr2['refmag'][0, :], randomn,
                                           calstr2['zmax'][0, :],
                                           calstr['etamin'][0], calstr['n0'][0],
                                           calstr2['volume'][0, :], calstr2['zrange'][0, :],
                                           calstr2['zbinsize'][0],
                                           zredstr, maxchi=20.0,
                                           ab_use=calstr2['afterburner_use'][0, :])

        # These match the IDL values
        testing.assert_almost_equal(rmfitter(calstr['cmax'][0, :]), 1.9331937798956758)

        p0_cval = np.zeros(calstr['nodes'][0, :].size) + 2.0
        testing.assert_almost_equal(rmfitter(p0_cval), 317.4524284321642)

        cvals = rmfitter.fit(p0_cval)

        # This does not match the IDL output, because this is doing a lot
        # better job minimizing the function, at least in this test.
        # I hope this is just because of the size of the testbed, which is
        # really way too small for something like this.
        testing.assert_almost_equal(cvals, np.array([2.61657263, 2.20376531, 1.00663991]))

        # Now we have to check the fitting with the afterburner

        biasvals = np.zeros(rmfitter._corrnodes.size)
        eratiovals = np.ones(rmfitter._corrnodes.size)
        biasvals, eratiovals = rmfitter.fit_bias_eratio(cvals, biasvals, eratiovals)

        cvals = rmfitter.fit(cvals, biaspars=biasvals, eratiopars=eratiovals, afterburner=True)

        testing.assert_almost_equal(cvals, np.array([3.39002141, 1.74421087, 0.93541002]))
        testing.assert_almost_equal(biasvals, np.array([0.00896487, -0.02456343, 0.02006761]))
        testing.assert_almost_equal(eratiovals, np.array([1.49999937, 1.01673233, 0.65657318]))

    def test_redmagic_calibrate(self):
        """
        """

        np.random.seed(12345)

        file_path = 'data_for_tests'
        conf_filename = 'testconfig_redmagic.yaml'
        config = Configuration(os.path.join(file_path, conf_filename))

        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestRedmapper-')
        config.outpath = self.test_dir

        testgals = GalaxyCatalog.from_fits_file(os.path.join('data_for_tests', 'redmagic_test', 'redmagic_test_input_gals.fit'))

        testgals.add_fields([('mag', 'f4', 5), ('mag_err', 'f4', 5)])

        redmagic_cal = RedmagicCalibrator(config)
        # We have to have do_run=False here because we don't have a real
        # galaxy training set with associated zreds!
        redmagic_cal.run(gals=testgals, do_run=False)

        # Read in the calibrated parameters

        self.assertTrue(os.path.isfile(config.redmagicfile))

        cal = fitsio.read(config.redmagicfile, ext=1)

        # Check that they are what we think they should be
        # (these checks are arbitrary, just to make sure nothing has changed)

        testing.assert_almost_equal(cal['cmax'][0, :], np.array([1.14588386, 3.89420298, -0.36792211]))
        testing.assert_almost_equal(cal['bias'][0, :], np.array([-0.09999912, -0.04537928, 0.01599778]))
        testing.assert_almost_equal(cal['eratio'][0, :], np.array([1.4999998, 1.48021495, 0.50000003]))

        pngs = glob.glob(os.path.join(self.test_dir, '*.png'))
        self.assertEqual(len(pngs), 3)

        # This is a hack of the volume limit mask to change from the one used for
        # calibration to the one used for the run (which uses a different footprint
        # because of reasons)

        config_regular = Configuration(os.path.join(file_path, 'testconfig.yaml'))
        maskfile = config_regular.maskfile

        config = Configuration(redmagic_cal.runfile)

        cal, hdr = fitsio.read(config.redmagicfile, ext=1, header=True)
        config.maskfile = maskfile
        os.remove(cal['vmaskfile'][0].decode().rstrip())
        mask = VolumeLimitMask(config, cal['etamin'], use_geometry=True)

        # Now test the running, using the output file which has valid galaxies/zreds
        run_redmagic = RunRedmagicTask(redmagic_cal.runfile)
        run_redmagic.run()

        # check that we have a redmagic catalog
        rmcatfile = config.redmapper_filename('redmagic_%s' % ('highdens'))
        self.assertTrue(os.path.isfile(rmcatfile))

        # And a random catalog
        rmrandfile = config.redmapper_filename('redmagic_%s_randoms' % ('highdens'))
        self.assertTrue(os.path.isfile(config.redmapper_filename('redmagic_%s_randoms' % ('highdens'))))

        # And check that the plot is there

        pngs = glob.glob(os.path.join(self.test_dir, '*.png'))
        self.assertEqual(len(pngs), 4)

        # And that we have the desired number of redmagic and randoms
        red_cat = GalaxyCatalog.from_fits_file(rmcatfile)
        rand_cat = GalaxyCatalog.from_fits_file(rmrandfile)

        self.assertEqual(rand_cat.size, red_cat.size * 10)

        # And confirm that all the randoms are in the footprint
        zmax = mask.calc_zmax(rand_cat.ra, rand_cat.dec)
        self.assertTrue(np.all(rand_cat.z < zmax))

    def setUp(self):
        self.test_dir = None

    def tearDown(self):
        if self.test_dir is not None:
            if os.path.exists(self.test_dir):
                shutil.rmtree(self.test_dir, True)


if __name__=='__main__':
    unittest.main()
