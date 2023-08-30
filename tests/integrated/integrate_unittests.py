# integrate_unittests.py
#
# Test Suite for the conintegrate function
# of the Picard Lefshetz module. 

import unittest
from plind.integrate.core import *
from plind.interpolate.core import *
import pl_testfunctions as plfun
import numpy as np

TOL= 10**-7

class TestConintegrate1D(unittest.TestCase):

	def test_trivial(self):
	# Integrates the zero function over the real line
		line= np.linspace(-1000, 1000, 10000)
		spln, spln_der, u_vals = spline1d(line)
		self.assertTrue(conintegrate(lambda x: (0+0*1j)*x, spln, spln_der, u_vals)< TOL)

	def test_line(self):
	# Integrates the function f(x)=x over [-1,1]
		line= np.linspace(-1, 1, 1000)
		spln, spln_der, u_vals = spline1d(line)
		self.assertTrue(conintegrate(plfun.real_line, spln, spln_der, u_vals)< TOL)

	def test_parabola(self):
	# Integrates the function f(x)=-x^2+3 over [-1,1]
		line= np.linspace(-1, 1, 1000)
		spln, spln_der, u_vals = spline1d(line)
		result = conintegrate(lambda x: -x**2+3, spln, spln_der, u_vals)
		self.assertTrue((result < 16/3+TOL) and (result > 16/3 -TOL))

	def test_gaussian(self):
	# Integrates the gaussian function over the real line
		line= np.linspace(-100, 100, 100)
		spln, spln_der, u_vals = spline1d(line)
		result = conintegrate(lambda x: np.exp(-x**2), spln, spln_der, u_vals)
		self.assertTrue((result < np.sqrt(np.pi)+TOL) and (result > np.sqrt(np.pi) -TOL))

	def test_complex(self):
	# Integrates x^2 as a complex function from [-i, i]
		line = np.linspace(-1,1, 10000)*1j
		spln, spln_der, u_vals = spline1d(line)
		result = conintegrate(lambda x: x**2, spln, spln_der, u_vals)
		real_correct = (result.real < TOL) and (result.real >  -TOL)
		imag_correct = (result.imag < -2/3+TOL) and (result.imag > -2/3-TOL)
		self.assertTrue(real_correct and imag_correct)

if __name__ == '__main__':
	unittest.main()

