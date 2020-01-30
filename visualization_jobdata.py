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
expfun = lambda x, lamb: 1j*lamb*(x**3/3-x)
start_time = 0.0
end_time = 100
Npts = 100
Nint = 100

#init_contour = np.exp(1j*0)*(1/np.tan(np.linspace(2*np.pi, -eps, Npts, endpoint=False)/2) + 0.0j)
#init_contour=np.exp(1j*np.linspace(0, 2*np.pi, Npts))
init_contour=np.linspace(-5,5, Npts)

# ------------------------------------------------------------------------------
# Detailed Contour Convergence for a given lambda.
#
# If the critical points of the function are known, this will also plot the
# Lefshetz thimbles so the progress towards them can be visualized.
# ------------------------------------------------------------------------------

# PARAMETERS TO DEFINE (!!!) #
lamb = 1
critpts = [-1, 1]
domain = [-10, 10]

# Generate model and perform descent
plind = plmodel(init_contour, expfun, expargs=[lamb])
plind.descend(start_time, end_time)

soln= plind.get_solution()
msg=soln.sget_message()

print(msg)

# -------------------------------------------------------------------------------
# Open up Job's data for the Airy function and plot
#
# -------------------------------------------------------------------------------
fpath="/home/etyhurst/picardlef/jobpl/Picard_Lefschetz_Integrator/code/data5/simplices.bin"
contourjob=np.fromfile(fpath, dtype=np.double)
print(np.size(contourjob))
contourjob=np.reshape(contourjob, [int(np.size(contourjob)/2),2])
contourjob=np.unique(contourjob, axis=0)
print(contourjob)
print(contourjob[:,0])

line = plind.get_contour()
trajectory = plind.get_trajectory()
hfun = plind.get_morse()
gradh = plind.get_grad()

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
im = ax.pcolormesh(u, v, hfun(z).real)
fig.colorbar(im, ax=ax)
for p0 in critpts:
    # Plot critical points
    ax.scatter(p0.real, p0.imag, color='w')
    # Plot contours of expfun.imag passing through p0, p1
    ax.contour(u, v, expfun(z,lamb).imag, [expfun(p0,lamb).imag], colors='w', alpha=0.5)
# Plot progress towards Lefschetz thimble
for k in np.linspace(0, trajectory.shape[0], trajectory.shape[0]//10, endpoint=False):
    k = int(k)
    ax.plot(trajectory[k].real[1:-1], trajectory[k].imag[1:-1], 'ro-')

# Plot final contour
dn = 5
# ax.quiver(u[::dn, ::dn], v[::dn, ::dn], gradh(z[::dn, ::dn], lamb).real, gradh(z[::dn, ::dn], lamb).imag, scale=140, color='k')
ax.plot(line.real[1:-1], line.imag[1:-1], 'bo-')
ax.plot(contourjob[:,0], contourjob[:,1], 'go-')
plt.show()
