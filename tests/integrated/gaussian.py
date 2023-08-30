from testfunction import *
import numpy as np


def Gauss_exp(z, ndim):
    gauss_sum=0
    if (ndim ==1):
        gauss_sum=z**2
    else: 
        for i in range(ndim):
            gauss_sum= gauss_sum + z[i]**2
    return 1j*(gauss_sum)

def Gauss_gradh(z, ndim):
    return 2*1j*np.conj(z)

def Gauss_answer(ndim):
    return (np.pi*1j)**(ndim/2)

def Gauss_thimble(x, ndim):
    return x+1j*x

def Gauss_thimble_dist(point, ndim):
    return np.linalg.norm(point.real - point.imag)

def Gaussian(ndim):
    """Returns the Gaussian integral of n-dimensions as a test function.
       
       Parameters
       ----------
        ndim: integer
            The dimension of the Gaussian.
            
       Returns
       ----------
        GaussianTestFn: plind.tests.testfunction Object
            A test function structure for testing plind.py

       See Also
       --------
           
        testfunction: A Test Function with known Lefshetz Thimbles, gradient, and integral over R^n.    
    """
    expfun = Gauss_exp
    gradh = Gauss_gradh
    expargs = [ndim]
    ndim = ndim
    thimbles = Gauss_thimble
    thimble_dist = Gauss_thimble_dist
    integral = Gauss_answer
        
    return testfunction(expfun, gradh, expargs, ndim, thimbles, thimble_dist, integral)
    