# Integrated tests for picard lefshetz

import sys
sys.path.append("..")

import unittest
from plind import *
import pl_testfunctions as plfun
import numpy as np

TOL= 10**-7

class integratedTestplmodel(unittest.TestCase):

   def test_gaussian_2D(self):
   # Integrates the gaussian function over the real line
      base_array = np.linspace(-100,100, 100)
      zero_array = np.zeros(np.shape(base_array))
      w,x,y,z = np.meshgrid(base_array, base_array, zero_array, zero_array)
      gauss_soln = sqrt((2*pi)**2)
      func = lambda x: (1/2)*(x[0]**2+x[1]**2)
      model = plmodel([

      spln, spln_der, u_vals = spline1d(line)
      self.assertTrue(conintegrate(lambda x: (0+0*1j)*x, spln, spln_der, u_vals)< TOL)

