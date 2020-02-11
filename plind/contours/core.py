import numpy as np

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
