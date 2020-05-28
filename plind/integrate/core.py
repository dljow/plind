#!/usr/bin/python3
# conintegrate.py - performs contour integral
import numpy as np
from scipy.interpolate import splprep, splev
from scipy.integrate import simps, quadrature, fixed_quad
from ..projection import *
import quadpy 

def conintegrate(f, contour, args=[], order=3):
    scheme = quadpy.nsimplex.grundmann_moeller(contour.ndim, order)
    simps = np.stack(contour.points[contour.simplices], axis=-2)
    val = scheme.integrate(lambda x:  f(x, *args), simps)

    # estimate error
    scheme = quadpy.nsimplex.grundmann_moeller(contour.ndim, order+1)
    order_up = scheme.integrate(lambda x: f(x, *args), simps)
    return np.sum(val), np.abs(np.sum(val)-np.sum(order_up))

# def _vol(points):
#     """Computes the volume of the (ndim)-simplex"""
#     ndim = points.shape[-1]
#     mat = (points - points[0])[1:]
#     # the sqrt(mat@mat) ensures that the orientation of all the simplices are consistent, but the overall integral can be off by an overall minus sign
#     # TODO: fix the overall minus sign possibly by tracking orientations.
#     if (np.shape(mat) == (1,1)) or (np.shape(mat) == (1,)):
#         dV = np.sqrt(np.squeeze(mat)**2)
#     else:
#         dV = np.sqrt(np.linalg.det(mat@mat))/np.math.factorial(ndim)
#     return dV
#
# def conintegrate(f, contour, args=[]):
#     integral = 0j
#     for simp in contour.simplices:
#         mid = np.sum(contour.points[simp], 0)/len(simp)
#         integral += f(mid, *args)*_vol(contour.points[simp])
#     return integral
