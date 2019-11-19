## Test module for the contour functions
#
# Tests the division of contour parts against a
# few simple cases.
import sys
sys.path.append("..")

import unittest
from core import *
import numpy as np


points = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
test_contour = contour(points=points)


class TestContour(unittest.TestCase):

  def test_init_contour(self):
  # Tests the setup of a set of points into a triangulation 
     test_contour.init_contour(points)
     
  # Tests that the shapes are as expected. 
     self.assertTrue(np.array_equal(points, test_contour.points))
     self.assertTrue(np.shape(test_contour.simplices), [3,2])
     self.assertTrue(np.shape(test_contour.edges), [2,5])

  def test_get_edgelengths(self):
  # Tests the edgelength function 
    test_contour.init_contour(points)

    lengths = test_contour.get_edgelengths()
    true_lengths= np.array([1,1,1,1,np.sqrt(2)])
    print(true_lengths)
    self.assertTrue(set(true_lengths)== set(lengths))

  def test_split_edges(self):
    self.assertTrue(False)

  def test_remove_points(self):
    self.assertTrue(False)


  def test_refine_edges(self):
    self.assertTrue(False)

if __name__ == '__main__':
        unittest.main()
