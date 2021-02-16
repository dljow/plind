import numpy as np


class testfunction:
    """A Test Function with known Lefshetz Thimbles, gradient, and integral over R^n. 
           
           Attributes
           ----------
           expfun: nfunction
            The function in the exponential of the oscillatory integrand. E.g. if the integrand
            is of the form exp(iS), then expfun = iS. expfun should take as its first argument a
            complex vector, z, in C^ndim (that is z is an ndarray<complex> of length ndim). It may take
            any number of additional arguments: expfun(z, *expargs).

           gradh: function
            The complex gradient of the Morse function, h = Real(iS). That is grad = dh/dz. The
            gradient is required to determine the direction in which to deform the contour. The
            arguments of grad should be the same as the arguments as expfun: grad(z, *expargs).
            
           expargs: array_like
             Additional arguments to be passed to expfun and grad. 
             
           ndim: int
             The dimension of the complex space, i.e. C^ndim.
           
           thimbles: function
             A function from R^ndim => C^ndim that parameterizes the thimbles. 
             
           thimble_dist: function
             A function from C^ndim => R that given a point purported to be on a Lefshetz thimble,
             returns an error on the point to the thimble. 
             
           integral: function
            The value of the integral for the given expargs. 
    """
    
    def __init__(self, expfun, gradh, expargs, ndim, thimbles, thimble_dist, integral):
        self.expfun = expfun
        self.gradh = gradh
        self.expargs = expargs
        if np.shape(self.expargs) == ():
            self.expargs = [self.expargs]

        self.ndim = ndim
        self.thimbles = thimbles
        self.thimble_dist = thimble_dist
        self.integral = integral 