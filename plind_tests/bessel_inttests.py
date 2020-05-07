## bessel_inttests.py
#
# Checks that the Picard-Lefshetz integrator
# can successfully integrate Bessel functions. 
#
# Can return analytics on the speed as well. 

import sys
sys.path.append("..")

import unittest
from plind.plmodel import *
import pl_testfunctions as plfun
import numpy as np
import time
import matplotlib.pyplot as plt

TOL= 10**-2

class TestBessel(unittest.TestCase):

        def test_j0(self):
            contour = np.linspace(-800, 1000, 500)
            j0 = lambda x:  1j*x-np.log(x)
            model =plmodel(contour, j0)
            model.descend(0, 0.5)
            model.integrate()
            ans=model.get_integral()
            print(ans)
            analytic=3.1416
            self.assertTrue(abs(ans-analytic) < TOL)

        def test_j0_powerseries(self):
            contour = np.linspace(-800, 1000, 500)
            j0 = lambda x:  1j*x-np.log(x)
            model =plmodel(contour, j0)
            model.descend(0, 0.5)
            model.integrate()
            ans=model.get_integral()
            print(ans)
            analytic=3.1416
            self.assertTrue(abs(ans-analytic) < TOL)



if __name__ == '__main__':
        unittest.main()


