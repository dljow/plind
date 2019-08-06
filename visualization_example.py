# visualization_example.py

"""Generic visualization script for a chosen example function."""

import numpy as np
from plind.plmodel import plmodel
import pl_testfunctions as plfun
import matplotlib.pyplot as plt

eps = np.finfo(float).eps

# ------------------------------------------------------------------------------
# Define the example function to integrate.
#
# Airy, Gaussian, and Witten's Airy are all defined in plfun.
# ------------------------------------------------------------------------------

# PARAMETERS TO DEFINE
expfun = plfun.AiryWittenexp
start_time = 0.0
end_time = 0.1
Npts = 100
Nint = 1000

init_contour = np.exp(1j*0)*(1/np.tan(np.linspace(2*np.pi, -eps, Npts, endpoint=False)/2) + 0.0j)

# ------------------------------------------------------------------------------
# Detailed Contour Convergence for a given lambda.
#
# If the critical points of the function are known, this will also plot the
# Lefshetz thimbles so the progress towards them can be visualized.
# ------------------------------------------------------------------------------

# PARAMETERS TO DEFINE (!!!) #
lamb = 1
critpts = [-1j, 1j]
domain = [-5, 5]

# Generate model and perform descent
plind = plmodel(init_contour, expfun, expargs=[lamb])
plind.descend(start_time, end_time)

line = plind.get_contour()
trajectory = plind.get_trajectory()
hfun = plind.get_morse()

# Plot final contour and trajectory
N = 100
u, v = np.meshgrid(np.linspace(domain[0], domain[1], N), np.linspace(domain[0], domain[1], N), indexing='ij')
z = u + 1j*v

fig, ax = plt.subplots()
ax.set_title('$h(z)$, $\\lambda$ = {}'.format(lamb))
ax.set_xlabel("Re($z$)")
ax.set_ylabel("Im($z$)")
ax.set_ylim(domain[0], domain[1])
ax.set_xlim(domain[0], domain[1])
# Plot h
im = ax.pcolormesh(u, v, hfun(z, lamb).real)
fig.colorbar(im, ax=ax)
for p0 in critpts:
    # Plot critical points
    ax.scatter(p0.real, p0.imag, color='w')
    # Plot contours of expfun.imag passing through p0, p1
    ax.contour(u, v, expfun(z, lamb).imag, [expfun(p0, lamb).imag], colors='w', alpha=0.5)
# Plot progress towards Lefschetz thimble
for k in np.arange(1, trajectory.shape[0], trajectory.shape[0]//5):
    ax.plot(trajectory[k].real[1:-1], trajectory[k].imag[1:-1], 'ro-')
# Plot final contour
ax.plot(line.real[1:-1], line.imag[1:-1], 'bo-')
plt.show()
