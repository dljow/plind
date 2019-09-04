#!/usr/bin/python3
# conintegrate.py - performs contour integral
import autograd.numpy as np
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

def conintegrate(f, contour_spline, contour_spline_der, spline_param, integrator=fixed_quad, Nint=200):
    integrand_R = lambda x: np.real( f(contour_spline(x)) * contour_spline_der(x) )
    integrand_I = lambda x: np.imag( f(contour_spline(x)) * contour_spline_der(x) )

    result_R = integrator(integrand_R, spline_param[0], spline_param[-1], n=Nint)
    result_I = integrator(integrand_I, spline_param[0], spline_param[-1], n=Nint)

    integral_R = result_R[0]
    error_R = result_R[1]

    integral_I = result_I[0]
    error_I = result_I[1]

    integral = integral_R + 1j*integral_I
    return integral
