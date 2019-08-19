#!/usr/bin/python3
import numpy as np
from scipy.interpolate import splprep, splev
from ..projection import *

def spline1d(line):
    pts = plane_to_sphere(line)
    tck, u = splprep(pts, s=0)

    line_map = lambda x: sphere_to_plane(splev(x, tck))
    line_tan = lambda x: sphere_to_plane_vec(splev(x, tck), splev(x, tck, der=1))
    return line_map, line_tan, u
