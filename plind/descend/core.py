import numpy as np
from ..projection import *

def identify_pole(points, gradh, dt, expargs, thresh=10):
    # performs Euler time step
    if np.any((np.abs(gradh(np.array(points), *expargs))/np.mean(np.abs(gradh(np.array(points), *expargs))))>thresh):
              print("A pole was found")

def flow(points, gradh, dt, expargs=[]):
    # performs Euler time step
   # identify_pole(points, gradh, dt, expargs)
    return np.array(points)+gradh(np.array(points), *expargs)*dt




# def flow_eq(time, line, gradh, args=[], anchor=[]):
#     """    Flow_eq defines the flow equation used to solve for the trajectory of the contour.
#     It returns the component of the gradient of the Morse function perpendicular to the
#     contour. It does this by taking in a contour in the complex plane, projecting to the
#     Riemann sphere, computing the gradient, and then projecting back into the complex
#     plane.
#
#     Parameters
#     ----------
#     time: float
#          the time at which the gradient is requested.
#
#     line: np.ndarray
#          the manifold to flow along the gradient.
#
#     gradh: func
#          the gradient of the Morse function h.
#
#     args: array-like (optional)
#          the arguments to the gradient of the Morse function, if needed.
#
#     Returns
#     -------
#     dydt: np.ndarray
#         the perpendicular gradient at all the points in line. """
#     n = line.shape[0]//2
#
#     u = line[:n]
#     v = line[n:2*n]
#     # Project onto Riemann sphere (embedded in 3D (x, y, z) coordinates)
#     eps = np.finfo(float).eps  # machine epsilon
#     denom = 1+u**2+v**2
#
#     x, y, z = plane_to_sphere(u+1j*v)
#
#     # Compute tangent vector in 3D
#     dx = (np.concatenate([x[1:], x[:1]]) - np.concatenate([x[-1:], x[:-1]]))/2
#     dy = (np.concatenate([y[1:], y[:1]]) - np.concatenate([y[-1:], y[:-1]]))/2
#     dz = (np.concatenate([z[1:], z[:1]]) - np.concatenate([z[-1:], z[:-1]]))/2
#
#     # Compute gradient in (u, v) space, then project onto Riemann sphere
#     # grad is the slowest thing right now
#     grad = gradh(u+1j*v, *args)
#
#     gradx, grady, gradz = plane_to_sphere_vec(u+1j*v, grad)
#
#     # Compute the perpendicular component of the gradient
#     mag = (gradx*dx + grady*dy + gradz*dz)/(dx**2+dy**2+dz**2)
#     gradperpx = gradx - mag*dx
#     gradperpy = grady - mag*dy
#     gradperpz = gradz - mag*dz
#
#
#     # Project back onto the complex plane
#     gradperp = sphere_to_plane_vec([x,y,z], [gradperpx,gradperpy,gradperpz])
#
#     fu = gradperp.real
#     fv = gradperp.imag
#     # set velocity of anchor to 0
#     fu[anchor] = 0.0
#     fv[anchor] = 0.0
#
#     # Change of variable ( t -> s/(1-s) ) to integrate to infinite time in finite parameter
#     fu *= (1-time)**(-2)
#     fv *= (1-time)**(-2)
#
#     dydt = np.concatenate((fu, fv))
#
#     return dydt
#
# def tot_speed(time, line, gradh, term_frac_eval=0.25, args=[]):
#     """Compute the total speed for the slowest `term_frac_eval` fraction of points."""
#     percentile = term_frac_eval*100
#     n = line.shape[0]//2
#
#     dydt = flow_eq(time, line, gradh, args)
#
#     dudt = dydt[:n]
#     dvdt = dydt[n:2*n]
#
#     slowest = (dudt**2+dvdt**2) < np.percentile((dudt**2+dvdt**2), percentile)
#
#     dudt = dudt[slowest]
#     dvdt = dvdt[slowest]
#
#     gradtot = np.sqrt(np.sum(dudt**2) + np.sum(dvdt**2))
#
#     return gradtot
#
# def terminal_cond(time, line, gradh, term_tol=np.inf, term_frac_eval=0.25, args=[]):
#     """Terminal condition event tracker for solve_ivp. Whenever gradtot is slower than a set value, return 0.
#     """
#     gradtot = tot_speed(time, line, gradh, term_frac_eval, args)
#     if gradtot <= term_tol:
#         return 0
#     else:
#         return 1
