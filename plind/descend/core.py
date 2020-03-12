import numpy as np
from ..projection import *

def identify_pole(points, gradh, dt, expargs, thresh=10):
    # performs Euler time step
    if np.any((np.abs(gradh(np.array(points), *expargs))/np.mean(np.abs(gradh(np.array(points), *expargs))))>thresh):
              print("A pole was found")

def _euler(points, gradh, dt, expargs):
    return points + dt * gradh(points, expargs)

def _rk4(points, gradh, dt, expargs):
    #
    k1 = dt * gradh(points, expargs)
    k2 = dt * gradh(points + k1/2, expargs)
    k3 = dt * gradh(points + k2/2, expargs)
    k4 = dt * gradh(points + k3, expargs)
    return points + (k1 + 2*k2 + 2*k3 + k4)/6

def flow(points, gradh, dt, expargs=[]):
    # performs Euler time step
   # identify_pole(points, gradh, dt, expargs)
    #flow_method = _euler
    flow_method = _rk4
    return flow_method(points, gradh, dt, *expargs)
