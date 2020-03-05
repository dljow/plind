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
