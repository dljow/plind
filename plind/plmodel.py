import autograd.numpy as np
from autograd import elementwise_grad as egrad
from scipy.integrate import solve_ivp, fixed_quad
from .plexception import *
from .integrate import conintegrate
from .poles import *
from .interpolate import *
from .descend import *
from .contour import *

DIVERGE = 10**9  # Divergence threshold. functions that evaluate to higher than this are considered poles.

class plmodel:
    """some documentation."""

    def __init__(self, contour, expfun, grad=None, expargs=[]):
        self.contour = contour
        self.expfun = expfun
        self.grad = grad
        self.expargs = expargs
        if np.shape(self.expargs) == ():
            self.expargs = [self.expargs]

        self.ndim = contour.simplices.shape[1]-1  # ndim = number of vertices in simplex minus 1
        self.trajectory = np.array([contour])
        self.integral = None
        self.critpts = []
        self.poles = []   # Identifies regions of the contour that may contain poles.

    # Simple functions for retrieving attributes
    def get_contour(self):
        return self.contour

    def get_expfun(self):
        return self.expfun

    def get_integral(self):
        return self.integral

    def get_critpts(self):
        return self.critpts

    # Functions for getting things that are derived from the attributes
    def get_trajectory(self):
        return self.trajectory

    # def get_contour_spline(self):
    #    return spline1d(self.contour.points)

    def get_poles(self):
        """Identifies parts of the domain that may be poles."""
        return self.poles, self.intfun(self.poles, *self.expargs)

    def get_intfun(self):
        """Return integrand function, i.e. np.exp(self.expfun(z, args=self.expargs))."""
        def intfun(z, *args):
            return np.exp(self.expfun(z, *args))
        return intfun


    def get_morse(self):
        """Return the morse function, i.e. the real part of expfun."""
        def morse(z, *args):
            return np.real(self.expfun(z, *args))
        return morse


    def get_grad(self, dx=1e-6):
        """Return self.grad. If self.grad is none, returns numerical gradient computed from self.expfun."""
        if self.grad is None:
            morse = self.get_morse()

            morse_grad = egrad(morse)

            def auto_grad(z, *args):
                # gradRe = np.real(morse_grad(z, *args))
                # gradIm = np.imag(morse_grad(z, *args))
                return -np.conj(morse_grad(z, *args))

            return auto_grad
        else:
            return self.grad

    # Functions for performing the PL integration
    def descend(self, dt, Nstep, delta):
        gradh = self.get_grad()

        i = 0
        while i < Nstep:
            # perform euler
            self.contour.points = flow(self.contour.points, gradh, dt, expargs=self.expargs)  # perform euler pushing

            # remove points
            bad_points = []  # find the points to remove
            self.contour.remove_points(bad_points)

            # refine mesh
            self.contour.refine_edges(delta)

            # add new contour to trajectory
            self.trajectory = np.append(self.trajectory, self.contour)

            i += 1

    # Function for integrating over contour
    def integrate(self, integrator=fixed_quad, Nint=200, intfun=None):
        if intfun is None:
            self.intfun = self.get_intfun()
        else:
            self.intfun = intfun

        self.integral = conintegrate(self.intfun, self.contour, args=self.expargs)

        # self.contour_spline, self.contour_spline_der, self.contour_spline_param = self.get_contour_spline()
        #
        # # Identify poles:
        # xvals = self.contour_spline(self.contour_spline_param)
        # eval = self.intfun(xvals, *self.expargs)
        #
        # if np.sum((abs(eval) > DIVERGE)) > 0:
        #         self.poles = xvals[(abs(eval) > DIVERGE)]
        #         raise PoleError("Poles were identified in the interpolated contour. Check self.poles.")
        #
        # intfun_wrapped = lambda z: self.intfun(z, *self.expargs)
        # self.integral = conintegrate(intfun_wrapped, self.contour_spline, self.contour_spline_der, self.contour_spline_param, integrator, Nint=Nint)
