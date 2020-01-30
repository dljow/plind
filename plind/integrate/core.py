#!/usr/bin/python3
# conintegrate.py - performs contour integral
import numpy as np
from scipy.interpolate import splprep, splev
from scipy.integrate import simps, quadrature, fixed_quad
from ..projection import *

def conintegrate(f, contour, args=[]):
    int = 0
    for simp in contour.simplices:
        mid = np.sum(contour.points[simp], 0)/len(simp)
        # compute volume element (this is very janky and not general at all)
        [p0, p1, p2] = contour.points[simp]
        base = np.sqrt(np.sum((p1-p0)**2))
        mid_pt = (p1+p0)/2
        height = np.sqrt(np.sum((p2-mid_pt)**2))
        #base = np.sum((p1-p0)*np.conj(p1-p0))
        #mid_pt = (p1+p0)/2
        #height = np.sum((p2-mid_pt)*np.conj(p2-mid_pt))
        vol = 0.5*base*height

        int += f(mid, *args)*vol
    return int

#def conintegrate(f, line, args=[], Nint=1000):
#    pts = p.plane_to_sphere(line)
#    tck, u = splprep(pts, s=0)  # it is not clear what s means, but s = None (default) breaks things
#    param_grid = np.linspace(u[0], u[-1], Nint)
#    vec = splev(param_grid, tck)
#    pts = vec[0]+1j*vec[1]
#    dpts = splev(param_grid, tck, der=1)
#    deriv = dpts[0]+1j*dpts[1]
#    integral = simps(f(pts, *args)*deriv, x=param_grid)
#    return integral
