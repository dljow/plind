import numpy as np

def flow_eq(time, line, gradh, args=[]):
    n = line.size//2

    u = line[:n]
    v = line[n:2*n]

    # Project onto Riemann sphere (embedded in 3D (x, y, z) coordinates)
    eps = np.finfo(float).eps  # machine epsilon
    denom = 1+u**2+v**2

    x = 2*u/denom
    y = 2*v/denom
    z = (u**2+v**2-1)/denom

    # Compute tangent vector in 3D
    dx = (np.concatenate([x[1:], x[:1]]) - np.concatenate([x[-1:], x[:-1]]))/2
    dy = (np.concatenate([y[1:], y[:1]]) - np.concatenate([y[-1:], y[:-1]]))/2
    dz = (np.concatenate([z[1:], z[:1]]) - np.concatenate([z[-1:], z[:-1]]))/2

    # Compute gradient in (u, v) space, then project onto Riemann sphere
    # (grad is the slowest thing right now)
    grad = gradh(u+1j*v, args)

    gradx = (2*(1-u**2+v**2)/denom**2 * grad.real) - (4*u*v/denom**2 * grad.imag)
    grady = (-4*u*v/denom**2 * grad.real) + (2*(1+u**2-v**2)/denom**2 * grad.imag)
    gradz = (4*u/denom**2 * grad.real) + (4*v/denom**2 * grad.imag)

    # Compute the perpendicular component of the gradient
    mag = (gradx*dx + grady*dy + gradz*dz)/(dx**2+dy**2+dz**2)
    gradperpx = gradx - mag*dx
    gradperpy = grady - mag*dy
    gradperpz = gradz - mag*dz

    # Project back onto the complex plane
    gradperpu = gradperpx / (1-z+eps) + gradperpz * x / (1-z+eps)**2
    gradperpv = gradperpy / (1-z+eps) + gradperpz * y / (1-z+eps)**2

    fu = gradperpu
    fv = gradperpv

    # Change of variable ( t -> s/(1-s) ) to integrate to infinite time in finite parameter
    fu *= (1-time)**(-2)
    fv *= (1-time)**(-2)

    return np.concatenate((fu, fv))
