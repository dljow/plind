import numpy as np
from scipy.misc import derivative
from .solution import solution
from .conintegrate import conintegrate

class plmodel:
    """some documentation."""

    def __init__(self, contour, expfun, grad=None, expargs=[]):
        self.contour = contour
        self.expfun = expfun
        self.grad = grad
        self.expargs = []

        self.solution = solution()
        self.integral = 0
        self.critpts = []


    def get_contour(self):
        return self.contour


    def get_trajectory(self):
        pass


    def get_expfun(self):
        return self.expfun


    def get_intfun(self):
        """Return integrand function, i.e. np.exp(self.expfun(z, args=self.expargs))."""
        def intfun(z, args=self.expargs):
            return np.exp(self.expfun(z, *args))
        return intfun


    def get_morse(self):
        """Return the morse function, i.e. the real part of expfun."""
        def morse(z, args=self.expargs):
            return self.expfun(z, *args).real
        return morse


    def get_grad(self, dx=1e-6):
        """Return self.grad. If self.grad is none, returns numerical gradient computed from self.expfun."""
        if self.grad is None:
            morse = self.get_morse()

            def num_grad(z):
                gradRe = -derivative(morse, x0=z, dx=dx)
                gradIm = derivative(morse, x0=z, dx=dx)
                return gradRe.real + 1j*gradIm.imag

            return num_grad
        else:
            return self.grad


    def get_integral(self):
        return self.integral


    def get_critpts(self):
        return self.critpts


    def descend(self):
        pass


    def integrate(self, Nint=1000):
        intfun = self.get_intfun()
        self.integral = conintegrate(intfun, self.contour, Nint=Nint)
