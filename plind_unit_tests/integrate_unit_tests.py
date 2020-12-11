## Test module for the integration functions
#
# Test that given a predictable contour to integrate
# over, the integration gives the correct answer. 
#
# Since quadpy is already tested, the real case need not
# be tested in an involved manner. The complex case is more
# pressing. 
#

import sys
sys.path.append("..")

import unittest
from plind.integrate import conintegrate
import numpy as np

NMAX = 5 # Test up to this dimension for generic ND test cases. 
TOL = 10**-7 # Numerical integration error 

class TestIntegrate(unittest.TestCase):
  def test_real_fn_2D(self):
    """Test the numerical integration of sin(x)+sin(y) 
       over [[-2pi, 2pi], [-2pi, 2pi]]""" 
    
    
    self.assertTrue(False)
    
  def test_real_fn_3D(self):
    self.assertTrue(False)
    
  def test_pure_complex_cnst_nD(self):
    self.assertTrue(False)
    
  def test_rotated_complex_cnst_nD(self):
    self.assertTrue(False)
    
  def test_pure_complex_polyn_nD(self):
    self.assertTrue(False)
    
  def test_rotated_complex_polyn_nD(self):
    self.assertTrue(False)
    
  def test_Gaussian_nD(self):
    self.assertTrue(False)
  
    