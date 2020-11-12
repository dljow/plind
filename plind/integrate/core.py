#!/usr/bin/python3
# conintegrate.py - performs contour integral
import numpy as np
from ..projection import *
from ..helpers import *
#import quadpy
from sympy import Rational as frac
from math import factorial as fact

# this bit of code is taken directly from quadpy. How does credit work here?
def get_vol(simplex):
    # Compute the volume via the Cayley-Menger determinant
    # <http://mathworld.wolfram.com/Cayley-MengerDeterminant.html>. One advantage is
    # that it can compute the volume of the simplex indenpendent of the dimension of the
    # space in which it is embedded.

    # compute all edge lengths
    edges = np.subtract(simplex[:, None], simplex[None, :])
    ei_dot_ej = np.einsum("...k,...k->...", edges, edges)

    j = simplex.shape[0] - 1
    a = np.empty((j + 2, j + 2) + ei_dot_ej.shape[2:], dtype=complex)
    a[1:, 1:] = ei_dot_ej
    a[0, 1:] = 1.0
    a[1:, 0] = 1.0
    a[0, 0] = 0.0

    a = np.moveaxis(a, (0, 1), (-2, -1))
    det = np.linalg.det(a)

    vol = np.sqrt((-1.0) ** (j + 1) / 2 ** j / fact(j) ** 2 * det)
    return vol


def grundmann_moeller_integrate(f, contour, order):
    s = order
    n = contour.ndim
    simps = np.stack(contour.points[contour.simplices], axis=-2)

    d = 2 * s + 1
    exponents = get_all_exponents(n + 1, s)
    data = [
        (
            frac(
                (-1) ** i * 2 ** (-2 * s) * (d + n - 2 * i) ** d,
                fact(i) * fact(d + n - i),
            ),
            np.array(
                [
                    [frac(2 * p + 1, d + n - 2 * i) for p in part]
                    for part in exponents[s - i]
                ]
            ),
        )
        for i in range(s + 1)
    ]
    points, weights = untangle(data)
    weights /= sum(weights)

    flt = np.vectorize(float)
    simplex = np.asarray(simps)
    x = np.dot(simplex.T, flt(points).T)
    vol = get_vol(simplex)

    fx = np.asarray(f(x))
    assert (
            x.shape[1:] == fx.shape[-len(x.shape[1:]):]
        ), "Illegal shape of f(x) (expected (..., {}), got {})".format(
            ", ".join([str(k) for k in x.shape[1:]]), fx.shape
        )
    return vol*np.dot(fx, flt(weights))

def conintegrate(f, contour, args=[], order=3):
    val = grundmann_moeller_integrate(lambda x:  f(x, *args), contour, order)
    order_up = grundmann_moeller_integrate(lambda x:  f(x, *args), contour, order+1)
    return np.sum(val), np.abs(np.sum(val)-np.sum(order_up))

# def conintegrate(f, contour, args=[], order=3):
#     scheme = quadpy.nsimplex.grundmann_moeller(contour.ndim, order)
#     simps = np.stack(contour.points[contour.simplices], axis=-2)
#     val = scheme.integrate(lambda x:  f(x, *args), simps)
#
#     # estimate error
#     scheme = quadpy.nsimplex.grundmann_moeller(contour.ndim, order+1)
#     order_up = scheme.integrate(lambda x: f(x, *args), simps)
#     return np.sum(val), np.abs(np.sum(val)-np.sum(order_up))

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
