# visualization_example.py

"""Generic visualization script for a chosen example function."""

import numpy as np
from plind.plmodel import plmodel
import pl_testfunctions as plfun
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from scipy.special import jv
from scipy.optimize import curve_fit

eps = np.finfo(float).eps

# ------------------------------------------------------------------------------
# Define the example function to integrate.
#
# Airy, Gaussian, and Witten's Airy are all defined in plfun.
# ------------------------------------------------------------------------------

def corner_fit(corners, xrange):
    """Fit a corner with a parabola"""
    
    parabola=lambda x, a,b,c: a*x**2+b*x+c
    popt, cov=curve_fit(parabola, np.real(corners), np.imag(corners))
    return xrange+parabola(xrange, *popt)*1j
    

# PARAMETERS TO DEFINE
expfun = lambda x, lamb: 2*lamb*np.sinh(x)+x
start_time = 0.0
end_time = 0.15
Npts = 300
Nint = 100

init_contour = np.exp(1j*0)*(1/np.tan(np.linspace(2*np.pi, -np.pi, Npts, endpoint=False)/2) + 0.0j)
#init_contour=np.exp(1j*np.linspace(0, 2*np.pi, Npts))

# Use the contour from 3.5 of Witten
cont_1=np.linspace(-0.3,-0.01,10) #np.exp(1j*0)*(1/np.tan(np.linspace(2*np.pi, -np.pi, Npts/2, endpoint=False)/2) + 0.0j)
cont_corner=corner_fit([-0.015,-0.05+0.005*1j, 1j*0.015], np.linspace(-0.015, 0,8))
print(cont_corner)
cont_1=np.append(cont_1, cont_corner)
cont_2=1j*np.linspace(cont_1[-1].imag+0.01,2*np.pi-cont_1[-1].imag-0.01, Npts)
cont_3= np.flip(cont_1)+1j*2*np.pi #np.flip(np.exp(1j*0)*(1/np.tan(np.linspace(-np.pi, -eps, Npts/2, endpoint=False)/2) + 0.0j))
init_contour=np.append(np.append(cont_1, cont_2), cont_3)
#print(init_contour)
#print(np.sum(np.isinf(init_contour)))
eval=np.exp(expfun(init_contour,1))
#print(eval)
#print(np.sum(eval>0))
#print(np.sum(np.isinf(eval)))


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
plind.integrate()
integral=plind.get_integral()
msg=soln.sget_message()
error=np.abs(integral/(2*np.pi*1j)-jv(2, 1))

print("Error in integral")
print(error)
print("Integral")
print(integral)
print("bessel")
print(jv(2,1))
print(msg)

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
im = ax.pcolormesh(u, v, hfun(z).real, norm=SymLogNorm(1))
fig.colorbar(im, ax=ax)
for p0 in critpts:
    # Plot critical points
    ax.scatter(p0.real, p0.imag, color='w')
    # Plot contours of expfun.imag passing through p0, p1
    ax.contour(u, v, expfun(z,lamb).imag, [expfun(p0,lamb).imag], colors='w', alpha=0.5)
# Plot progress towards Lefschetz thimble
for k in np.linspace(0, trajectory.shape[0], trajectory.shape[0]//50, endpoint=False):
    k = int(k)
    ax.plot(trajectory[k].real[1:-1], trajectory[k].imag[1:-1], 'ro-')

# Plot final contour
dn = 5
# ax.quiver(u[::dn, ::dn], v[::dn, ::dn], gradh(z[::dn, ::dn], lamb).real, gradh(z[::dn, ::dn], lamb).imag, scale=140, color='k')
ax.plot(line.real[1:-1], line.imag[1:-1], 'bo-')
plt.show()

# Compare integration to the actual Spherical Bessel function in question


