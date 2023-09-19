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

import pytest
import numpy as np
import contour from plind

@pytest.fixture
def example_contour():
    points = np.array([[0, 0], [0, 1], [1, 0], [1, 1], [1,-1]])
    return contour(points=points)

def test_init_contour():
    return