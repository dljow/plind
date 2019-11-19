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
    self.assertTrue(set(true_lengths)== set(lengths))

  def test_split_edges(self):
    # Tests if bad edges are appropriately split
    bad_edges =[[0,3]] # diagonal edge in block

    test_contour.init_contour(points)
  
    test_contour.split_edges(bad_edges)
    new_points = test_contour.points.tolist()
    new_edges = test_contour.edges.tolist()
    new_simplices= test_contour.simplices.tolist()

    # Check the appropriate components have been added
    self.assertTrue([0.5,0.5] in new_points)
    self.assertTrue([1,4] in new_edges)
    self.assertTrue([2,4] in new_edges)
    self.assertTrue([0,1,4] in new_simplices)
    self.assertTrue([0,2,4] in new_simplices)
    self.assertTrue([1,3,4] in new_simplices)

    # Check the appropriate components have been removed
    self.assertTrue([0,3] not in new_edges)
    self.assertTrue([0,1,3] not in new_simplices)
    self.assertTrue([0,2,3] not in new_simplices)

  def test_remove_points(self):
  # Tests removal of points
    bad_point= 0
    test_contour.init_contour(points)
    test_contour.remove_points([bad_point])
    new_points = test_contour.points.tolist()
    new_edges = test_contour.edges.tolist()
    new_simplices= test_contour.simplices.tolist()

    # point is actually removed
    self.assertTrue(bad_point not in new_points)

    # edges are back indexed appropriately
    self.assertTrue([0,2] in new_edges)
    self.assertTrue([1,2] in new_edges)
    self.assertTrue(len(new_simplices)==0)


  def test_refine_edges(self):
  # Tests overall refinement of an edge, this is an integrated
  # test of all of the above. 
    delta = np.sqrt(2)-0.001 

    test_contour.init_contour(points)
    test_contour.refine_edges(delta)

    new_points = test_contour.points.tolist()
    new_edges = test_contour.edges.tolist()
    new_simplices= test_contour.simplices.tolist()

    # The diagonal edge is no longer present
    self.assertTrue([0,3] not in new_edges)

    # The too-large simplices are no longer present
    self.assertTrue([0,1,3] not in new_simplices)
    self.assertTrue([0,2,3] not in new_simplices)

    # The midpoint is in the array 
    self.assertTrue([0.5, 0.5] in new_points) 

    # The new edges have been added
    self.assertTrue([3,4] in new_edges)
    self.assertTrue([1,3] in new_edges)
    self.assertTrue([0,4] in new_edges)
    self.assertTrue([2,4] in new_edges)

    # The new simplices have been added
    self.assertTrue([0,1,4] in new_simplices)
    self.assertTrue([0,2,4] in new_simplices)
    self.assertTrue([1,3,4] in new_simplices)
    self.assertTrue([2,3,4] in new_simplices)

if __name__ == '__main__':
        unittest.main()
