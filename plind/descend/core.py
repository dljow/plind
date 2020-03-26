import numpy as np
from ..projection import *

def identify_pole(points, gradh, dt, expargs, thresh=10):
    # performs Euler time step
    if np.any((np.abs(gradh(np.array(points), *expargs))/np.mean(np.abs(gradh(np.array(points), *expargs))))>thresh):
              print("A pole was found")

def _euler(points, gradh, dt, expargs):
    return points + dt * gradh(points, expargs), dt

def _rk4(points, gradh, dt, expargs):
    #
    k1 = dt * gradh(points, expargs)
    k2 = dt * gradh(points + k1/2, expargs)
    k3 = dt * gradh(points + k2/2, expargs)
    k4 = dt * gradh(points + k3, expargs)
    return points + (k1 + 2*k2 + 2*k3 + k4)/6, dt

def _rkf45(points, gradh, dt, expargs):
    #HARDCODED STEP-WISE TOLERANCE!!
    tol = 1e-5

    k1 = dt * gradh(  points, expargs)
    k2 = dt * gradh(  points + k1 * 1/4, expargs)
    k3 = dt * gradh(  points + k1 * 3/32      + k2 * 9/32, expargs)
    k4 = dt * gradh(  points + k1 * 1932/2197 + k2 * -7200/2197 + k3 * 7296/2197, expargs)
    k5 = dt * gradh(  points + k1 * 439/216   + k2 * -8         + k3 * 3680/513   + k4 * -845/4104, expargs)
    k6 = dt * gradh(  points + k1 * -8/27     + k2 * 2          + k3 * -3544/2565 + k4 * 1859/4104   + k5 * -11/40, expargs)

    points_next     = points + k1 * 25/216    +                   k3 * 1408/2565  + k4 * 2197/4104   + k5 * -1/5
    points_next_alt = points + k1 * 16/135    +                   k3 * 6656/12825 + k4 * 28561/56430 + k5 * -9/50 + k6 * 2/55

    R = np.abs(points_next - points_next_alt) / dt
    delta = np.min((np.divide(tol,2*R))**(1/4))

    if np.any(R > tol):
        dt = delta*dt
        _rkf45(points, gradh, dt, expargs)
    return points_next, delta*dt

def _dop853(points, gradh, dt, expargs):
    k1 = dt * gradh(points, expargs)
    k2 = dt * gradh(points + k1/5, expargs)
    k3 = dt * gradh(points + k1*3/40 + k2*9/40, expargs)
    k4 = dt * gradh(points + k1*44/45 + k2*-56/15 + k3*32/9, expargs)
    k5 = dt * gradh(points + k1*19372/6561 + k2*-25360/2187 + k3*64448/6561 + k4*-212/729, expargs)
    k6 = dt * gradh(points + k1*9017/3168 + k2*-355/33 + k3*46732/5247 + k4*49/176 + k5*-5103/18656, expargs)
    k7 = dt * gradh(points + k1*35/384 + k3*500/1113 + k4*125/192 + k5*-2187/6784 + k6*11/84, expargs)
    points_next_1 = points + k7
    points_next_2 = points + k1*5179/57600 + k3*7571/16695 + k4*393/640 + k5*-92097/339200 + k6*187/2100 + k7/40

def flow(points, gradh, dt, expargs=[]):
    # performs Euler time step
   # identify_pole(points, gradh, dt, expargs)
    #flow_method = _euler
    flow_method = _rkf45
    return flow_method(points, gradh, dt, *expargs)
