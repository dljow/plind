# pt_population_inttests.py
#
# Test Suite for verifying the plind program
# does not overpopulate points at infinity
# and underpopulate the critical points

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


class TestInftyPts(unittest.TestCase):
        """Tests for ensuring the descend function does not have points that fly off to infinity"""

        def test_domainconverge(self):
        # Tests that the trajectory of the gradient flow remains roughly constant after an increase to 10000 iterates
           domain= np.linspace(-10,10, 100)
           domain_model = PLModel(domain, expfun, gradh, expargs=[lamb])
           step_range= np.linspace(minstep, maxstep, 10)
           domain_model.descend(0, minstep-100)
           last_contour= domain_model.get_contour()


           for end_time in step_range:
                   # check differences of the two endpoints
                   domain_model.contour = domain
                   domain_model.descend(0, end_time)
                   current_contour= domain_model.contour
                   diff = abs(current_contour[0] - last_contour[0])+ abs(current_contour[-1]-last_contour[-1])
                   self.assertTrue(diff < DIST)
                   last_contour= current_contour


        def test_critpt_pop(self):
        # Tests that the contour remains well-populated amoung critical points
           domain= np.linspace(-10,10, 100)
           domain_model = PLModel(domain, expfun, gradh, expargs=[lamb])
           domain_model.descend(0, minstep)

           line = domain_model.contour

           # keep count of total number of points near critical point:
           nearpoints= 0

           for critpt in critpts:
                    distances = np.abs(line-critpt)
                    countpts= np.sum(distances < CRITRAD)
                    nearpoints+=countpts

           # count the percentage of points near a critical point

           perc = 100*(nearpoints/np.size(line))
           print(perc)
           self.assertTrue(perc > PERCENT )



if __name__ == '__main__':
	unittest.main()
