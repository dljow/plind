from .testfunction import testfunction
import numpy as np


def Airy_exp(z, ndim):
    Airy_sum=0
    for i in range(ndim):
        Airy_sum= Airy_sum + z[i]**2
    return 1j*(Airy_sum)

def Airy_gradh(z, ndim):
    return 2*1j*np.conj(z)

def Airy_answer(ndim):
    return (np.pi*1j)**(ndim/2)

def Airy_thimble(x, ndim):
    return x+1j*x

def Airy_thimble_dist(point, ndim):
    return np.abs(point.real - point.imag)

def Airy(ndim):
    """Returns the Airy integral of n-dimensions as a test function.
       
       Parameters
       ----------
        ndim: integer
            The dimension of the Airy function.
            
       Returns
       ----------
        AiryTestFn: plind.tests.testfunction Object
            A test function structure for testing plind.py

       See Also
       --------
           
        testfunction: A Test Function with known Lefshetz Thimbles, gradient, and integral over R^n.    
    """
    expfun = Airy_exp
    gradh = Airy_gradh
    expargs = [ndim]
    ndim = ndim
    thimbles = Airy_thimble
    thimble_dist = Airy_thimble_dist
    integral = Airy_answer
        
    return testfunction(expfun, gradh, expargs, ndim, thimbles, thimble_dist, integral)
    
