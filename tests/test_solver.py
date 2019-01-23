from __future__ import division, absolute_import, print_function
from past.builtins import xrange

import unittest
import numpy.testing as testing
import numpy as np
import fitsio

import redmapper

class SolverNFWTestCase(unittest.TestCase):
    """
    Tests for redmapper zero-finding solver.
    """
    def runTest(self):
        """
        Run tests on redmapper.Solver
        """
        file_name = 'test_solver_data.fit'
        file_path = 'data_for_tests'

        data=fitsio.read('%s/%s' % (file_path,file_name),ext=1)

        #need to transpose cpars
        data[0]['CPARS'] = data[0]['CPARS'][::-1]

        # check some common errors...
        testing.assert_raises(ValueError,redmapper.Solver,data[0]['R0'],data[0]['BETA'],data[0]['UCOUNTS'][0:10],data[0]['BCOUNTS'],data[0]['R'],data[0]['W'],cpars=data[0]['CPARS'],rsig=data[0]['RSIG'])
        testing.assert_raises(ValueError,redmapper.Solver,data[0]['R0'],data[0]['BETA'],data[0]['UCOUNTS'],data[0]['BCOUNTS'][0:10],data[0]['R'],data[0]['W'],cpars=data[0]['CPARS'],rsig=data[0]['RSIG'])
        testing.assert_raises(ValueError,redmapper.Solver,data[0]['R0'],data[0]['BETA'],data[0]['UCOUNTS'],data[0]['BCOUNTS'],data[0]['R'][0:10],data[0]['W'],cpars=data[0]['CPARS'],rsig=data[0]['RSIG'])
        testing.assert_raises(ValueError,redmapper.Solver,data[0]['R0'],data[0]['BETA'],data[0]['UCOUNTS'],data[0]['BCOUNTS'],data[0]['R'],data[0]['W'][0:10],cpars=data[0]['CPARS'],rsig=data[0]['RSIG'])
        testing.assert_raises(ValueError,redmapper.Solver,data[0]['R0'],data[0]['BETA'],data[0]['UCOUNTS'],data[0]['BCOUNTS'],data[0]['R'],data[0]['W'],cpars=data[0]['CPARS'][0:1],rsig=data[0]['RSIG'])

        # and test the results
        solver=redmapper.Solver(data[0]['R0'],data[0]['BETA'],data[0]['UCOUNTS'],data[0]['BCOUNTS'],data[0]['R'],data[0]['W'],cpars=data[0]['CPARS'],rsig=data[0]['RSIG'])

        """
        solve_nfw() spits out:
        lambda,
        p_mem,
        wt = p_mem*theta^L*theta^R
        """
        lam,p,wt,rlambda,theta_r=solver.solve_nfw()

        testing.assert_almost_equal(lam,data[0]['LAMBDA'])
        testing.assert_array_almost_equal(p,data[0]['PVALS'])
        testing.assert_array_almost_equal(wt,data[0]['WTVALS'])
        testing.assert_almost_equal(rlambda,data[0]['R0']*(data[0]['LAMBDA']/100.)**data[0]['BETA'])
        testing.assert_almost_equal(theta_r,data[0]['THETA_R'],6)


if __name__=='__main__':
    unittest.main()
