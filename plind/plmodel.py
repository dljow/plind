import numpy as np
from .solution import solution

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

    def get_trajector(self):
        pass

    def get_expfun(self):
        return self.expfun

    def get_morse(self):
        """Return the morse function, i.e. the real part of expfun."""
        def morse(z, args=self.expargs):
            return self.expfun(z, *args).real
        return morse

    def get_grad(self):
        """Return self.grad. If self.grad is none, returns numerical gradient computed from self.expfun."""
        if self.grad is None:
            pass
        else:
            return self.grad

    def get_integral(self):
        return self.integral

    def get_critpts(self):
        return self.critpts

    def descend(self):
        pass

    def integrate(self):
        pass
