import numpy as np
import matplotlib.pyplot as plt

from ..plmodel import plmodel

def visualize_descent(plmodel=plmodel, domain = [-5,5], lognorm = False):
    morsefun = plmodel.get_morse()
    expfun = plmodel.expfun
    trajectory = plmodel.get_trajectory()

    init_contour = trajectory[0]
    line = plmodel.contour

    expargs = plmodel.expargs
    critpts = plmodel.critpts

    N = 100
    U = np.linspace(domain[0], domain[1], N)
    V = np.linspace(domain[0], domain[1], N)
    U, V = np.meshgrid(U, V, indexing='ij')
    Z = U + 1j*V

    fig, ax = plt.subplots()
    ax.set_title('$h(z)$, $\\args$ = {}'.format(expargs))
    ax.set_xlabel("Re($z$)")
    ax.set_ylabel("Im($z$)")
    ax.set_ylim(domain[0], domain[1])
    ax.set_xlim(domain[0], domain[1])
    # Plot h
    im = ax.pcolormesh(U, V, morsefun(Z, expargs))
    fig.colorbar(im, ax=ax)
    for p0 in critpts:
        # Plot critical points
        ax.scatter(p0.real, p0.imag, color='w')
        # Plot contours of expfun.imag passing through p0, p1
        ax.contour(U, V, expfun(Z, lamb).imag, [expfun(p0, expargs).imag], colors='w', alpha=0.5)
    # Plot progress towards Lefschetz thimble
    for k in np.linspace(0, trajectory.shape[0], trajectory.shape[0]//10, endpoint=False):
        k = int(k)
        ax.plot(trajectory[k].real[1:-1], trajectory[k].imag[1:-1], 'ro-')

    # Plot final contour
    dn = 5
    # ax.quiver(u[::dn, ::dn], v[::dn, ::dn], gradh(z[::dn, ::dn], lamb).real, gradh(z[::dn, ::dn], lamb).imag, scale=140, color='k')
    ax.plot(line.real[1:-1], line.imag[1:-1], 'bo-')
    plt.show()
