# integrated_test_fixed_params.py
#
# A generic integrated

import sys
sys.path.append("..")

import unittest
from plind.plmodel import *
from plind.contour_dict import equilateral_real, realcontour_nd, realcontour_1D
from plind.plexception.plexception import *
import gaussian as gauss

import numpy as np


def setup_descend(testfunction, contour, parameters):
    """Sets up the Picard-Lefshetz gradient flow
        and runs for the given parameters. This setup is
        run at the start of every test, that is every test is
        with a fresh run of the gradient flow.
        
        Parameters
        -------
            testfunction: plind.tests.testfunction Object
                 A test function structure for testing plind.py
                 
            parameters: dict
                 A dictionary containing the parameters passed
                 to descend()

        
        Returns
        -------
            plind: plind.plmodel Object
                 A plmodel object with the gradient flow desended
                 by the timesteps provided by dict. 
                
        See Also
        --------
           
        testfunction: A Test Function with known Lefshetz Thimbles, gradient, and integral over R^n.
        
        """    
    
    plind = plmodel(contour, testfunction.expfun, grad=testfunction.gradh, expargs=testfunction.expargs)
    plind.descend(parameters["delta"], parameters["thresh"], parameters["tmax"], parameters["dt_init"])
    return plind


class FixedParamsPlindTest(unittest.TestCase):
    THIMBLE_TOL=10**-2
    DELTA = 0.6
    NSTEP = 150
    DT_INIT = 1e-2
    THRESH =-7
    TMAX = DT_INIT*150
    
    def test_Gauss_1D_finds_contours(self):
        gauss_fn = gauss.Gaussian(1)
        param_dict = {
        "delta": self.DELTA,
        "nstep": self.NSTEP,
        "dt_init": self.DT_INIT,
        "thresh": self.THRESH,
        "tmax": self.TMAX
        }
        contour = realcontour_1D(10, (-1.5,1.5))
        
        plind = setup_descend(gauss_fn, contour, param_dict)
        
        contour = plind.contour
        points= contour.points
        
        for point in points:
            self.assertTrue(gauss_fn.thimble_dist(point,1)< self.THIMBLE_TOL)
        
if __name__ == '__main__':
	unittest.main()
        
        