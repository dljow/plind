#!/usr/bin/python3
# conintegrate.py - performs contour integral
import numpy as np
from ..helpers import *
from sympy import Rational as frac
from math import factorial as fact


def grundmann_moeller_integrate(f, contour, order):
    '''Implementation of the Grundmann-Moeller integration scheme from

        Source code: quadpy
        Author: Nico Schlomer
        Code version: 0.16.6
        Availability: https://github.com/nschloe/quadpy'''

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
    vol = contour.get_vols(imag=True)

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
