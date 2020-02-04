import numpy as np
import plind.contour as ctr
import itertools

# returns an equilateral triangle tesselation of R2
def equilateral_real(N, domain):
    contour = ctr.contour()
    grid_x, grid_y = np.meshgrid(np.linspace(domain[0], domain[1], N), np.linspace(domain[2], domain[3], N), indexing='xy')
    delta = abs(domain[1]-domain[0])/N
    offset = np.ones(grid_x.shape)*delta
    for i in np.arange(0, N):
        if (i+1) % 2 != 0:
            offset[i, :] = 0
    grid_x = grid_x + offset
    contour.points = np.array([grid_x.flatten(), grid_y.flatten()]).T
    edges = []
    for i in np.arange(0, N-1):
        if i % 2 == 0:
            edges.append([i*N, (i+1)*N])
            for j in np.arange(1, N):
                edges.append([i*N+j, i*N+j-1])
                edges.append([i*N+j, (i+1)*N+j-1])
                edges.append([i*N+j, (i+1)*N+j])
        elif i % 2 != 0:
            edges.append([i*N+N-1, (i+1)*N+N-1])
            for j in np.arange(0, N-1):
                edges.append([i*N+j, (i+1)*N+j])
                edges.append([i*N+j, i*N+j+1])
                edges.append([i*N+j, (i+1)*N+j+1])
    for j in np.arange(0, N-1):
        edges.append([(N-1)*N+j, (N-1)*N+j+1])
    simplices = []
    for i in np.arange(0, N-1):
        for j in np.arange(0, N-1):
            if i % 2 == 0:
                simplices.append([i*N+j, i*N+j+1, (i+1)*N+j])
                simplices.append([(i+1)*N+j, (i+1)*N+j+1, i*N+j+1])
            elif i % 2 != 0:
                simplices.append([i*N+j, (i+1)*N+j, (i+1)*N+j+1])
                simplices.append([i*N+j, i*N+j+1, (i+1)*N+j+1])
    contour.edges = np.array(edges)
    contour.simplices = np.array(simplices)
    contour.ndim = contour.simplices.shape[1]-1
    return contour



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
