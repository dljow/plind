# integrated_test_fixed_params.py
#
# A generic integrated

import unittest
from plind.plmodel import *
from plind.contour_dict import equilateral_real, real_contour_nd, real_contour_1d
from plind.plexception.plexception import *
import gaussian as gauss
import itertools

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
            plind: plind.PLModel Object
                 A PLModel object with the gradient flow desended
                 by the timesteps provided by dict. 
                
        See Also
        --------
           
        testfunction: A Test Function with known Lefshetz Thimbles, gradient, and integral over R^n.
        
        """    
    
    plind = PLModel(contour, testfunction.expfun, grad=testfunction.gradh, expargs=testfunction.expargs)
    plind.descend(parameters["delta"], parameters["thresh"], parameters["tmax"], parameters["dt_init"])
    return plind

def setup_integrate(testfunction, contour, parameters):
    plind = setup_descend(testfunction, contour, parameters)
    plind.integrate()
    return plind


class FixedParamsPlindTest(unittest.TestCase):
    

## TOLERANCE PARAMETERS ##
    THIMBLE_TOL=10**-2
    CONSTANT_TOL = 5*10**-2
    INTEGRAL_TOL = 5*10**-3


    DELTA = 0.6
    NSTEP = 140
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
        contour = real_contour_1d(10, (-1.5,1.5))
        
        plind = setup_descend(gauss_fn, contour, param_dict)
        
        contour = plind.contour
        points= contour.points
        
        for point in points:
            self.assertTrue(gauss_fn.thimble_dist(point,1)< self.THIMBLE_TOL)
            
    def test_Gauss_1D_constant_eval(self):
        gauss_fn = gauss.Gaussian(1)
        param_dict = {
        "delta": self.DELTA,
        "nstep": self.NSTEP,
        "dt_init": self.DT_INIT,
        "thresh": self.THRESH,
        "tmax": self.TMAX
        }
        contour = real_contour_1d(10, (-1.5,1.5))
        
        plind = setup_descend(gauss_fn, contour, param_dict)
        
        contour = plind.contour
        points= contour.points
        
        for pair in itertools.combinations(points, 2):
            pairdiff = np.abs(gauss_fn.expfun(pair[0],1).imag - gauss_fn.expfun(pair[1],1).imag)
            self.assertTrue(pairdiff < self.CONSTANT_TOL)
            
    def test_Gauss_1D_integral(self):
        gauss_fn = gauss.Gaussian(1)
        param_dict = {
        "delta": self.DELTA,
        "nstep": self.NSTEP,
        "dt_init": self.DT_INIT,
        "thresh": self.THRESH,
        "tmax": self.TMAX
        }
        contour = real_contour_1d(10, (-1.5,1.5))
        
        plind = setup_integrate(gauss_fn, contour, param_dict)
        #print(plind.integral[1])
        self.assertTrue(np.abs(plind.integral[0]- gauss_fn.integral(1)) < self.INTEGRAL_TOL)
        #self.assertTrue(np.abs(plind.integral[0]- gauss_fn.integral(1)) < plind.integral[1])
        
            
    def test_Gauss_nD_finds_contours(self):
        param_dict = {
        "delta": self.DELTA,
        "nstep": self.NSTEP,
        "dt_init": self.DT_INIT,
        "thresh": self.THRESH,
        "tmax": self.TMAX
        }
        
        for n in [2,3]:
            gauss_fn = gauss.Gaussian(n)
            domain = tuple(np.ndarray.flatten(np.transpose(np.reshape(np.repeat((-1.5, 1.5), n), [2,n]))))
            contour = real_contour_nd(10, domain)
        
            plind = setup_descend(gauss_fn, contour, param_dict)
        
            contour = plind.contour
            points= contour.points
        
            for point in points:
                self.assertTrue(gauss_fn.thimble_dist(point,n)< self.THIMBLE_TOL)
        
if __name__ == '__main__':
	unittest.main()
        
        