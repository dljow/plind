## Test module for the contour functions
#
# Tests the contour function against a much more
# involved case. This is a key integration test
# for all the contour functions
#
# This contour test is applied to the following shape:
#
#  1 - - - e1 - - - 3
#  -              - -
#  -     s0     -   -
#  -          -     -
#  e0      e4       e2
#  -     -          -
#  -   -     s1     - 
#  - -              - 
#  0 - - - e3 - - - 2
#   -               -
#    -     s2       -
#     -             -
#      e5           e6
#         -         -
#           -       -
#             -     - 
#               -   4
#                   
#
# The points are indexed as I've drawn them.
# I make no guarantee of the edges/simplices 
# being properly indexed

import unittest
from plind.contour import *
import numpy as np

points = np.array([[0, 0], [0, 1], [1, 0], [1, 1], [1,-1]])
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
  # Tests that the edge length function 
  # properly returns the simplex edge lengths. 
    test_contour.init_contour(points)

    lengths = test_contour.get_edgelengths()
    lengths = np.ndarray.flatten(lengths)
    true_lengths= [1.,1.,1.,1.,1.,1.,np.sqrt(2), np.sqrt(2),np.sqrt(2)] # shared edge will appear twice
                                                       # since list is given by simplex
    self.assertTrue(sorted(true_lengths)== sorted(lengths.tolist()))

  def test_remove_points(self):
    test_contour.init_contour(points)
    
  # Tests removal of points
    bad_point_index= [4]
    bad_point =[1,-1]

    test_contour.init_contour(points)
    test_contour.remove_points([bad_point_index])
    new_points = test_contour.points.tolist()
    new_edges = test_contour.edges.tolist()
    new_simplices= test_contour.simplices.tolist()

    # point is actually removed
    self.assertTrue(bad_point not in new_points)

    # edges are back indexed appropriately
    self.assertTrue([[4, 2], [4, 0], [2, 0]] not in new_edges)
    self.assertTrue([[3, 1], [3, 0], [1, 0]] in new_edges)
    self.assertTrue([1,2] not in new_edges)
    self.assertTrue(len(new_simplices)==2)


  def test_refine_delta(self):
  # Tests overall refinement of an edge, this is an integrated
  # test of all of the above. 
    delta = np.sqrt(2)-0.001 

    # Exclude the last point for simplicity
    test_contour.init_contour(points[:-1])
    test_contour.refine_edges(delta)

    new_points = test_contour.points.tolist()
    new_edges = test_contour.edges.tolist()
    new_edges = np.reshape(new_edges, [int(np.size(new_edges)/2),2]).tolist()
    new_simplices= test_contour.simplices.tolist()

    # The diagonal edge is no longer present
    self.assertTrue([0,3] not in new_edges)

    # The too-large simplices are no longer present
    self.assertTrue([0,1,3] not in new_simplices)
    self.assertTrue([0,2,3] not in new_simplices)

    # The midpoint is in the array 
    self.assertTrue([0.5, 0.5] in new_points) 

    # The new edges have been added
    self.assertTrue([4,3] in new_edges)
    self.assertTrue([3,1] in new_edges)
    self.assertTrue([4,0] in new_edges)
    self.assertTrue([4,2] in new_edges)

    # The new simplices have been added
    self.assertTrue([4,1,0] in new_simplices)
    self.assertTrue([4,2,0] in new_simplices)
    self.assertTrue([4,3,1] in new_simplices)
    self.assertTrue([4,2,3] in new_simplices)

if __name__ == '__main__':
        unittest.main()
