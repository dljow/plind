#!/usr/bin/python3
# conintegrate.py - performs contour integral
import numpy as np
from scipy.interpolate import splprep, splev
from scipy.integrate import simps, quadrature, fixed_quad
from ..projection import *

def _vol(points):
    """Computes the volume of the (ndim)-simplex"""
    ndim = points.shape[-1]
    mat = (points - points[0])[1:]
    # the sqrt(mat@mat) ensures that the orientation of all the simplices are consistent, but the overall integral can be off by an overall minus sign
    # TODO: fix the overall minus sign possibly by tracking orientations.
    if (np.shape(mat) == (1,1)) or (np.shape(mat) == (1,)):
        dV = np.sqrt(np.squeeze(mat)**2)
    else:
        dV = np.sqrt(np.linalg.det(mat@mat))/np.math.factorial(ndim)
    return dV

def conintegrate(f, contour, args=[]):
    integral = 0j
    for simp in contour.simplices:
        mid = np.sum(contour.points[simp], 0)/len(simp)
        integral += f(mid, *args)*_vol(contour.points[simp])
    return integral

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

# def conintegrate(f, contour_spline, contour_spline_der, spline_param, integrator=fixed_quad, Nint=400):
#     """ Integrates the function f over the manifold line.
#     Parameters
#     ----------
#     f: function
#        the function to be integrated.
#
#     contour_spline: np.ndarray
#          the manifold to integrate the function over.
#
#     contour_spline_der: np.ndarray
#          the derivative of the contour spline
#
#     spline_param: array-like
#          the arguments to Morse function, if needed.
#
#     integrator: function (optional)
#          the integrator to us
#
#     Nint: integer (optional)
#          number of points to integrate over
#
#     Returns
#     -------
#     dydt: np.ndarray.3
#
#         the perpendicular gradient at all the points in line. """
#
#
#     integrand_R = lambda x: np.real( f(contour_spline(x)) * contour_spline_der(x) )
#     integrand_I = lambda x: np.imag( f(contour_spline(x)) * contour_spline_der(x) )
#
#     result_R = integrator(integrand_R, spline_param[0], spline_param[-1], n=Nint)
#     result_I = integrator(integrand_I, spline_param[0], spline_param[-1], n=Nint)
#
#     integral_R = result_R[0]
#     error_R = result_R[1]
#
#     integral_I = result_I[0]
#     error_I = result_I[1]
#
#     integral = integral_R + 1j*integral_I
#     return integral
