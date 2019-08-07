import numpy as np
from scipy.misc import derivative
from scipy.integrate import solve_ivp

from .integrate import conintegrate
from .descend import *

from .solution import solution

class plmodel:
    """some documentation."""

    def __init__(self, contour, expfun, grad=None, expargs=[]):
        self.contour = contour
        self.expfun = expfun
        self.grad = grad
        self.expargs = expargs
        if np.shape(self.expargs) == ():
            self.expargs = [self.expargs]

        self.solution = solution()
        self.integral = None
        self.critpts = []

        self._dimsize = np.shape(self.contour)

    # Simple functions for retrieving attributes
    def get_contour(self):
        return self.contour

    def get_expfun(self):
        return self.expfun

    def get_solution(self):
        return self.solution

    def get_integral(self):
        return self.integral

    def get_critpts(self):
        return self.critpts

    # Functions for getting things that are derived from the attributes
    def get_trajectory(self):
        return self.trajectory


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
    def descend(self, start_time, end_time, term_frac_eval = 0.25, term_percent = 0.1):
        gradh = self.get_grad()
        y0 = np.concatenate((self.contour.real, self.contour.imag))
        init_speed = tot_speed(start_time, y0, gradh, term_frac_eval, self.expargs)
        term_tol = init_speed*term_percent
        flow = lambda t, y: flow_eq(t, y, gradh, self.expargs)
        term_cond = lambda t, y: terminal_cond(t, y, gradh, term_tol, term_frac_eval, self.expargs)
        term_cond.terminal = True

        self.solution = solve_ivp(fun=flow, t_span=(start_time, end_time), y0=y0, method = 'BDF', vectorized = 'True')

        # time steps
        self._t = self.solution.t # time steps

        # stopping time, if exists
        self._t_events = self.solution.t_events

        # trajectory in returned form
        self._y = self.solution.y
        self.trajectory = self._y[:self._dimsize[0]] + 1j*self._y[self._dimsize[0]:2*self._dimsize[0]]

        #get the last contour of the trajectory
        self.contour = self.trajectory[...,-1]

    def integrate(self, Nint=1000):
        self.intfun = self.get_intfun()
        self.integral = conintegrate(lambda z: self.intfun(z, *self.expargs), self.contour, Nint=Nint)
