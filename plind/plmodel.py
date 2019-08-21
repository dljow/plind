import numpy as np
from scipy.misc import derivative
from scipy.integrate import solve_ivp, simps, quad, quadrature, fixed_quad
from .plexception import * 
from .integrate import conintegrate
from .poles import *
from .interpolate import *
from .descend import *
from .solution import solution

DIVERGE= 10**9 # Divergence threshold. functions that evaluate to higher than this are considered poles. 

class plmodel:
    """some documentation."""

    def __init__(self, contour, expfun, grad=None, expargs=[]):
        self.contour = contour
        self.expfun = expfun
        self.grad = grad
        self.expargs = expargs
        if np.shape(self.expargs) == ():
            self.expargs = [self.expargs]

        self.solution = None  # Initialize to none, remind user to descend
        self.integral = None
        self.critpts = []
        self.poles= [] # Identifies regions of the contour that may contain poles. 

    # Simple functions for retrieving attributes
    def get_contour(self):
        return self.contour

    def get_expfun(self):
        return self.expfun

    def get_solution(self):
        if self.solution is None:
            raise DescendError("You must call .descend to get a solution")
        return self.solution

    def get_integral(self):
        return self.integral

    def get_critpts(self):
        return self.critpts

    # Functions for getting things that are derived from the attributes
    def get_trajectory(self):
        if solution is None:
            raise DescendError("You must call .descend to get a solution")
            return None
        else:
            return self.solution.get_trajectory()

    def get_contour_spline(self):
        return spline1d(self.contour)

    def get_poles(self):
        """Identifies parts of the domain that may be poles."""
        return self.poles, self.intfun(poles, *self.expargs)

    def get_intfun(self):
        """Return integrand function, i.e. np.exp(self.expfun(z, args=self.expargs))."""
        def intfun(z, *args):
            return np.exp(self.expfun(z, *args))
        return intfun


    def get_morse(self):
        """Return the morse function, i.e. the real part of expfun."""
        def morse(z, *args):
            return self.expfun(z, *args).real
        return morse


    def get_grad(self, dx=1e-6):
        """Return self.grad. If self.grad is none, returns numerical gradient computed from self.expfun."""
        if self.grad is None:
            morse = self.get_morse()

            def num_grad(z, *args):
                gradRe = -derivative(lambda z: morse(z, *args), x0=z, dx=dx)
                gradIm = derivative(lambda z: morse(z, *args), x0=z, dx=dx*1j)
                return gradRe.real + 1j*gradIm.imag

            return num_grad
        else:
            return self.grad

    # Functions for performing the PL integration
    def descend(self, start_time, end_time, term_frac_eval=0.25, term_percent=0.1):
        gradh = self.get_grad()
        y0 = np.concatenate((self.contour.real, self.contour.imag))
        init_speed = tot_speed(start_time, y0, gradh, term_frac_eval, self.expargs)
        term_tol = init_speed*term_percent

        def flow(t, y):
            return flow_eq(t, y, gradh, self.expargs)

        def term_cond(t, y):
            return terminal_cond(t, y, gradh, term_tol, term_frac_eval, self.expargs)
        term_cond.terminal = True

        self.solution = solution(solve_ivp(fun=flow, t_span=(start_time, end_time), y0=y0, method='BDF', vectorized='True', events=term_cond))
        self.contour = self.solution.get_contour()

    def integrate(self, integrator=fixed_quad, Nint=200):
        self.intfun = self.get_intfun()
        self.contour_spline, self.contour_spline_der, self.contour_spline_param = self.get_contour_spline()

        # Identify poles:
        xvals= self.contour_spline(self.contour_spline_param)
        eval= self.intfun(xvals, *self.expargs)

        if np.sum((abs(eval)> DIVERGE))> 0:
                self.poles=xvals[(abs(eval)> DIVERGE)]
                raise PoleError("Poles were identified in the interpolated contour. Check self.poles.")

        intfun_wrapped = lambda z: self.intfun(z, *self.expargs)
        self.integral = conintegrate(intfun_wrapped, self.contour_spline, self.contour_spline_der, self.contour_spline_param, integrator, Nint=Nint)


