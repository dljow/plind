import numpy as np
from ..contour import *
from scipy.spatial import Delaunay
import itertools

# returns contour for a Delaunay triangulation of Rn

def real_contour_nd(N, domain, rotation=0.):
    """Initialize a contour object over R^ndim based on the domain given, with N points along
       each edge.

           Parameters
           ----------

           N: integer
              The number of points along each edge, that is the domain runs from

           domain: (2*ndim) tuple of floats
              Describes the start, end of every dimension of the square, that is
              (x1min, x1max, x2min, x2max, ...)

           rotation: float in (0, 2pi)
              Rotates the contour into the complex plane.

           Returns
           -------
           plind.contour Object


           See Also
           --------

            __init__: Initializes a plind.contour given points, edges, and simplices.

    """
    # Construct the grid
    ndim = int(len(domain)//2)
    linspaces = []
    for i in np.arange(ndim):
        linspaces.append(np.linspace(domain[i], domain[i+1], N))
    R = np.meshgrid(*linspaces)

    # Flatten for use with Delaunay
    flattened_comps = []
    for comp in R:
        flattened_comps.append(comp.flatten())
    points = np.dstack(flattened_comps)[0]

    # Use Delaunay to get simplices
    tri = Delaunay(points)

    # Construct contour object
    ctr = contour(points*np.exp(1j*rotation)+0*1j, tri.simplices)
    return ctr



# returns regular grid for R1
def real_contour_1d(N, domain):
    points = np.linspace(domain[0], domain[1], N)
    simplices = np.array([[j, j+1] for j in np.arange(0, len(points)-1)])
    points = np.reshape(points, [N, 1])

    return contour(points, simplices)


def _rotate(contour, pivot, angle):
    assert np.isreal(angle), "Angle to rotate should be real."
    return (np.exp(1j*angle) * (contour - pivot)) + pivot

def compact_lens_1d(x0, Npts=51, Npts_vert=4, inner_domain=[-3, 3]):
    """Returns a sensible initial contour for expfun of type i*((x-x0)**2 / 2 + psi(x)), where psi(x) vanishes at +- infty"""
    assert np.isreal(x0), "Source position should be real."

    # generate the real projective line, removing the endpts because they blow up, and leaving space for x0 to be a pt
    line = np.linspace(-np.pi,np.pi,Npts+2-1)[1:-1]/2
    line = np.tan(line)

    # x0 should be between the ends
    assert x0 > line[0] and x0 < line[-1], "Source position should be inside the initial line. Increase Npts?"

    # insert x0
    line = np.insert(line, np.searchsorted(line,x0), x0)

    # rotate the contour around the pivot x0, leave inner_domain on the real line, and adding vertical contours
    if x0 <= inner_domain[0]:
        inner_domain[0] = x0
        contour_l = _rotate(line[line <= inner_domain[0]], x0, np.pi/4)
        contour_r = _rotate(line[line  > inner_domain[1]], x0, np.pi/4)
        contour_h = line[(line > inner_domain[0]) & (line <= inner_domain[1])]
        contour_v = np.linspace(contour_h[-1], contour_r[0], Npts_vert+2, endpoint=True)[1:-1]
        contour = np.concatenate((contour_l, contour_h, contour_v, contour_r))
    elif x0 > inner_domain[1]:
        inner_domain[1] = x0
        contour_l = _rotate(line[line <= inner_domain[0]], x0, np.pi/4)
        contour_r = _rotate(line[line  > inner_domain[1]], x0, np.pi/4)
        contour_h = line[(line > inner_domain[0]) & (line <= inner_domain[1])]
        contour_v = np.linspace(contour_l[-1], contour_h[0], Npts_vert+2, endpoint=True)[1:-1]
        contour = np.concatenate((contour_l, contour_v, contour_h, contour_r))
    elif (x0 > inner_domain[0]) & (x0 <= inner_domain[1]):
        contour_l = _rotate(line[line <= inner_domain[0]], x0, np.pi/4)
        contour_r = _rotate(line[line  > inner_domain[1]], x0, np.pi/4)
        contour_h = line[(line > inner_domain[0]) & (line <= inner_domain[1])]
        contour_vl = np.linspace(contour_l[-1], contour_h[0], Npts_vert//2+2, endpoint=True)[1:-1]
        if   Npts_vert%2 == 0:
            contour_vr = np.linspace(contour_h[-1], contour_r[0], Npts_vert//2+2, endpoint=True)[1:-1]
        elif Npts_vert%2 == 1:
            contour_vr = np.linspace(contour_h[-1], contour_r[0], Npts_vert//2+1+2, endpoint=True)[1:-1]
        contour = np.concatenate((contour_l, contour_vl, contour_h, contour_vr, contour_r))

    assert np.size(contour) == Npts + Npts_vert, "Total number of points not equal to Npts + Npts_vert."
    return contour
