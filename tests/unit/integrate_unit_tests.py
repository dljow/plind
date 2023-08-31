## Test module for the integration functions
#
# Test that given a predictable contour to integrate
# over, the integration gives the correct answer. 
#
# Since quadpy is already tested, the real case need not
# be tested in an involved manner. The complex case is more
# pressing. 
#

import unittest
from plind.integrate import conintegrate
import numpy as np
from plind.contour_dict import real_contour_nd, real_contour_1d

NMAX = 5 # Test up to this dimension for generic ND test cases. 
TOL2D = 10**-7 # Numerical integration error 
TOL3D = 10**-4 # Recommended to set this lower
TOLND = 10**-4

class TestIntegrate(unittest.TestCase):
  def test_real_fn_2D(self):
    """Test the numerical integration of sin(x)+sin(y) 
       over [[-2pi, 2pi], [-2pi, 2pi]]""" 
    
    contour = real_contour_nd(30, (-2*np.pi,2*np.pi,-2*np.pi,2*np.pi))
    integral = conintegrate(lambda xvec: np.sin(xvec[0])+np.sin(xvec[1]), contour) 
    self.assertTrue(integral[0] < TOL2D)
    
  def test_real_fn_3D(self):
    """Test the numerical integration of sin(x)+sin(y)+cos(z) 
     over [[-2pi, 2pi], [-2pi, 2pi], [-2pi, 2pi]]""" 
    contour = real_contour_nd(15, (-2*np.pi,2*np.pi,-2*np.pi,2*np.pi, -2*np.pi,2*np.pi))
    integral = conintegrate(lambda xvec: np.sin(xvec[0])+np.sin(xvec[1])+np.cos(xvec[2]), contour) 
    self.assertTrue(integral[0] < TOL3D)
    
  def test_pure_complex_cnst_nD(self):
    """Test the numerical integration of 1 
     over an N-dimensional square to NMAX 
     in pure complex space""" 
    self.assertTrue(False)
    
  def test_rotated_complex_cnst_nD(self):
    """Test the numerical integration of 1 
     over an N-dimensional square to NMAX, 
     rotated in the plane 45 degrees"""
        
    self.assertTrue(False)
    
  def test_pure_complex_polyn_nD(self):
    """Test the numerical integration of the polynomial
     3x**2 -1 over over an N-dimensional square to NMAX, 
     in pure complex space"""
    self.assertTrue(False)
    
  def test_rotated_complex_polyn_nD(self):
    """Test the numerical integration of the polynomial
     3x**2 -1 over over an N-dimensional square to NMAX, 
     rotated in the plane 45 degrees"""
    self.assertTrue(False)
    
  def test_Gaussian_nD(self):
    """Test integrating a Gaussian over its
     Lefshetz thimbles (known) to NMAX dimensions."""
    self.assertTrue(False)
  
 
if __name__ == '__main__':
        unittest.main()   