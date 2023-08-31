# Integrated tests for picard lefschetz

import unittest
from plind.plmodel import *
from plind.contour.core import *
import numpy as np
from plind.contour_dict import equilateral_real


TOL= 10**-7
dt=10**-2
Nstep= 180
delta= 1
distTOL= dt*delta*2

class integratedTestplmodel(unittest.TestCase):

   def test_gaussian_contour(self):
   # Integrates the gaussian function over the real line
      base_array = np.linspace(-1,1, 10)+0*1j
      w,x = np.meshgrid(base_array, base_array)
      w=w.flatten()
      x=x.flatten()
      gauss_soln = np.sqrt((2*np.pi)**2)
      gradh = lambda x,lamb: 2*1j*lamb*np.conj(x)
      func = lambda x,lamb: lamb*np.dot(x,x)
      #cont = contour()
      #cont.init_contour(np.transpose([w,x]))
      #print(cont.points)
      #print(cont.edges)
      cont = equilateral_real(50, (-1,1,-1,1))
      model = PLModel(cont, func, gradh, expargs=[1])
      model.descend(dt, Nstep, delta,-4)
      traj = model.trajectory
      final_cont= model.contour.points
      print(np.shape(final_cont))
      for point in final_cont:
          self.assertTrue(abs(point[0].real-point[0].imag) < distTOL)
          self.assertTrue(abs(point[1].real-point[1].imag) < distTOL)



   def test_gaussian_2D_integral(self):
   # Integrates the gaussian function over the real line
      base_array = np.linspace(-1,1, 20)+0*1j
      w,x = np.meshgrid(base_array, base_array)
      w=w.flatten()
      x=x.flatten()
      gauss_soln = np.pi*1j
      gradh = lambda x: np.transpose(np.array([-x[:,0].imag+ 1j*x[:,0].real, -x[:,1].imag +1j*x[:,1].real]))
      func = lambda x: (1/2)*(x[0]**2+x[1]**2)
      cont = contour()
      cont.init_contour(np.transpose([w,x]))
      model = PLModel(cont, func, gradh)
      model.descend(dt, Nstep, delta,-5)
      model.integrate()

      ans = model.integral
      print(ans)
      print(gauss_soln)

      self.assertTrue(abs(ans-gauss_soln)< TOL)



if __name__ == '__main__':
        unittest.main()
