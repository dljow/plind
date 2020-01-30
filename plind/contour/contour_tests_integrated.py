## Test module for the contour functions
#
# Tests the contour function against a much more
# involved case. This is a key integration test
# for all the contour functions
#
import sys
sys.path.append("..")

import unittest
from core import *
import numpy as np


points = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
test_contour = contour(points=points)


class TestContour(unittest.TestCase):