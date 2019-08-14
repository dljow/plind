#!/usr/bin/python3
# conintegrate.py - performs contour integral
import numpy as np
from scipy.interpolate import splprep, splev
from scipy.integrate import simps, quadrature, fixed_quad
from ..projection import *

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

def conintegrate(f, line, args=[], Nint=1000):
    pts = plane_to_sphere(line)
    tck, u = splprep(pts, s=0)

    line_map = lambda x: sphere_to_plane(splev(x, tck))
    line_tan = lambda x: sphere_to_plane_vec(splev(x, tck), splev(x, tck, der=1))

    integrand_R = lambda x: ( f(line_map(x), *args) * line_tan(x) ).real
    integrand_I = lambda x: ( f(line_map(x), *args) * line_tan(x) ).imag

    result_R = fixed_quad(integrand_R, u[0], u[-1], n=200)
    result_I = fixed_quad(integrand_I, u[0], u[-1], n=200)

    integral_R = result_R[0]
    error_R = result_R[1]

    integral_I = result_I[0]
    error_I = result_I[1]

    integral = integral_R + 1j*integral_I
    return integral
