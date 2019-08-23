# pole_inttests.py
#
# Integrated test Suite for verifying the process of 
# excising poles from integrals does not prevent
# the integral from giving the correct answer. 

import sys
sys.path.append("..")

import unittest
from plind.plmodel import *
import pl_testfunctions as plfun
import numpy as np
import time
import matplotlib.pyplot as plt

# ----------------------------------------------
# Parameters
#
# ----------------------------------------------

TOL= 10**-5 # tolerance demanded from previous functions
DIST= 10**-2 # sum total distance of the Lefshetz thimble endpoints from each other between the steps 

CRITRAD = 3 # radius around which we want the critical points to be populated
PERCENT = 50 # expected percent threshold of points within CRITRAD of the critical points 

# Maxiumum Step:
minstep= 1200 #the step size after which we expect convergence
maxstep= 2200 #the range endpoint (going past the convergence step)

# Function to integrate (DEFINE):
expfun= plfun.Airyexp
intfun= plfun.Airyint
exactfun= plfun.Airyexact
lamb = -0.2

critpts= [-np.sqrt(abs(lamb)), np.sqrt(abs(lamb))]

# No need to define these, they should be taken care of
hfun= lambda x,lamb: plfun.h(expfun, x, lamb)
gradh= lambda x,lamb: plfun.gradh(hfun, x, lamb)

# ----------------------------------------------
# Test Suite
#
# ----------------------------------------------  


class TestPoleFlagging(unittest.TestCase):

        def test_flag_pole_given_contour(self):
	# Tests that given a contour to integrate over specifically containing a
	# pole to make sure it is flagged. 
            expfun = lambda x: 1j*(1/(1+x))
            contour = np.linspace(-10, 10, 20)
            model = plmodel(contour, expfun)
            try:
                model.integrate()
                ans=model.get_integral()
            except PoleError:
                poles, vals= model.get_poles()
                self.assertTrue(-1 in poles)

if __name__ == '__main__':
        unittest.main()


