#!/usr/bin/python3
# conintegrate.py - performs contour integral
import numpy as np
from scipy.interpolate import splprep, splev
from scipy.integrate import simps

def conintegrate(f, line, Nint, args=[]):
    tck, u = splprep([line.real, line.imag], s=0)  # it is not clear what s means, but s = None (default) breaks things
    param_grid = np.linspace(u[0], u[-1], Nint)
    vec = splev(param_grid, tck)
    pts = vec[0]+1j*vec[1]
    dpts = splev(param_grid, tck, der=1)
    deriv = dpts[0]+1j*dpts[1]
    integral = simps(f(pts, *args)*deriv, x=param_grid)
    return integral
