# Integrated tests for picard lefshetz

import sys
sys.path.append("..")

import unittest
from plind import *
import numpy as np

TOL= 10**-7
dt=10**-6
Nstep=500
delta= 50

class integratedTestplmodel(unittest.TestCase):

   def test_gaussian_2D(self):
   # Integrates the gaussian function over the real line
      base_array = np.linspace(-100,100, 100)+0*1j
      w,x = np.meshgrid(base_array, base_array)
      gauss_soln = sqrt((2*pi)**2)
      gradh = lambda x: [2*x[0].real, -2*x[0].imag, 2*x[1].real, 2*x[1].imag]
      func = lambda x: (1/2)*(x[0]**2+x[1]**2)
      model = plmodel([np.flatten(w), np.flatten(x)], func, gradh=gradh)
      model.descend(dt, Nstep, delta)
      traj = model.trajectory
      print(traj)
      model.integrate()

      ans = model.integral
      print(ans)
      print(gauss_soln)

      self.assertTrue(abs(ans-gauss_soln)< TOL)

