#!/usr/bin/python3
# projection.py - projects from complex plane to riemann sphere and vice-versa.
import numpy as np

eps = np.finfo(float).eps

def plane_to_sphere(z):
    """
    Project a point from the complex plane to the Riemann sphere, embedded in R^3
    """
    u = z.real
    v = z.imag
    denom = 1+u**2+v**2

    X = 2*u / denom
    Y = 2*v / denom
    Z = (u**2 + v**2 - 1) / denom

    pt = np.array([X, Y, Z])

    return pt

def plane_to_sphere_vec(z, vec_z):
    """
    Pushforward a tangent vector from the complex plane to the Riemann sphere, embedded in R^3
    """
    u = z.real
    v = z.imag

    vec_u = vec_z.real
    vec_v = vec_z.imag

    denom = (1 + u**2 + v**2)**2

    vec_X = 2 * (1-u**2+v**2)/denom * vec_u - 4*u*v/denom * vec_v
    vec_Y = -4*u*v/denom * vec_u + 2 * (1+u**2-v**2)/denom * vec_v
    vec_Z = 4*u/denom * vec_u + 4*v/denom * vec_v

    vec = np.array([vec_X,vec_Y,vec_Z])
    return vec

def sphere_to_plane(pt):
    """
    Project a point from the Riemann sphere to the complex plane.
    """
    X, Y, Z = pt[0], pt[1], pt[2]
    u = X / (1-Z)
    v = Y / (1-Z)

    z = u + 1j*v

    return z

def sphere_to_plane_vec(pt, vec):
    """
    Pushforward a tangent vector from the Riemann sphere to the complex plane.
    """
    X, Y, Z = pt[0], pt[1], pt[2]
    vec_X, vec_Y, vec_Z = vec[0], vec[1], vec[2]

    vec_u = vec_X/(1-Z+eps) + vec_Z * X/(1-Z+eps)**2
    vec_v = vec_Y/(1-Z+eps) + vec_Z * Y/(1-Z+eps)**2

    vec_z = vec_u + 1j * vec_v
    return vec_z

def plane_to_cyl(z):
    """
    Project a point from the complex plane to a cylinder, wrapping around the real axis, embedded in R^3
    """
    u = z.real
    v = z.imag
    denom = 1+u**2

    X = 2*u / denom
    Y = v
    Z = (u**2 + - 1) / denom

    pt = np.array([X, Y, Z])

    return pt

def plane_to_cyl_vec(z, vec_z):
    """
    Pushforward a tangent vector from the complex plane to a cylinder, embedded in R^3
    """
    u = z.real
    v = z.real

    vec_u = vec_z.real
    vec_v = vec_z.imag

    denom = (1 + u**2)**2

    vec_X = 2 * (1-u**2)/denom * vec_u
    vec_Y = vec_v
    vec_Z = 4*u/denom * vec_u

    vec = np.array([vec_X,vec_Y,vec_Z])
    return vec

def cyl_to_plane(pt):
    """
    Project a point from the Riemann sphere to the complex plane.
    """
    X, Y, Z = pt[0], pt[1], pt[2]
    u = X / (1-Z)
    v = Y

    z = u + 1j*v

    return z

def cyl_to_plane_vec(pt, vec):
    """
    Pushforward a tangent vector from the Riemann sphere to the complex plane.
    """
    X, Y, Z = pt[0], pt[1], pt[2]
    vec_X, vec_Y, vec_Z = vec[0], vec[1], vec[2]

    vec_u = vec_X/(1-Z) + vec_Z * X/(1-Z)**2
    vec_v = vec_Y

    vec_z = vec_u + 1j * vec_v
    return vec_z
