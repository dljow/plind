import numpy as np
from scipy.misc import derivative
from scipy.integrate import solve_ivp

from .integrate import conintegrate
from .descend import flow_eq

from .solution import solution
from .conintegrate import conintegrate
from .descend.flow_equation import flow_eq


class plmodel:
    """some documentation."""

    def __init__(self, contour, expfun, grad=None, expargs=[]):
        self.contour = contour
        self.expfun = expfun
        self.grad = grad
        self.expargs = expargs

        self.solution = solution()
        self.integral = None
        self.critpts = []

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
        pass


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
    def descend(self, start_time, end_time):
        gradh = self.get_grad()
        self.solution = solve_ivp(fun=lambda t, y: flow_eq(t, y, gradh, self.expargs), t_span=(start_time, end_time), y0=np.concatenate((self.contour.real, self.contour.imag)))
        # self.contour = ? do something to get contour from solution


    def integrate(self, Nint=1000):
        intfun = self.get_intfun()
        self.integral = conintegrate(lambda z: intfun(z, *self.expargs), self.contour, Nint=Nint)
