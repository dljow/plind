from .testfunction import testfunction
import numpy as np

# Test function that is a variant on the Fresnel
# Integral. Only valid up to dimension 2. See
# https://arxiv.org/pdf/1201.1975.pdf (may extend
# to dimension 4) 


def Fresnel_exp(z, ndim):
    if ndim== 1:
        return 1j*z**2
    if ndim== 2:
        return 1j*(z[0]**2 - z[1]**2)

def Fresnel_gradh(z, ndim):
    if ndim== 1:
        return 2.*1j*np.conj(z)
    if ndim== 2:
        print(z[5])
        gradfres=np.transpose(np.array([2.*np.conj(z[:,0]), -2.*np.conj(z[:,1])]))
        print(gradfres[5])
        return gradfres
    
def Fresnel_answer(ndim):
    if ndim== 1:
        return np.sqrt(np.pi/2)
    if ndim== 2:
        return np.pi/2

def Fresnel_thimble(x, ndim):
    return None

def Fresnel_thimble_dist(point, ndim):
    return None

def Fresnel(ndim):
    """Returns the Fresnel integral of 1,2-dimensions as a test function.
       
       Parameters
       ----------
        ndim: integer
            The dimension of the Fresnel-like integral (1,2).
            
       Returns
       ----------
        FresnelTestFn: plind.tests.testfunction Object
            A test function structure for testing plind.py

       See Also
       --------
           
        testfunction: A Test Function with known Lefshetz Thimbles, gradient, and integral over R^n.    
    """
    expfun = Fresnel_exp
    gradh = Fresnel_gradh
    expargs = [ndim]
    ndim = ndim
    thimbles = Fresnel_thimble
    thimble_dist = Fresnel_thimble_dist
    integral = Fresnel_answer
        
    return testfunction(expfun, gradh, expargs, ndim, thimbles, thimble_dist, integral)
    
