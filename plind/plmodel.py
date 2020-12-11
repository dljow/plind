
import numpy as np
from scipy.integrate import fixed_quad
from copy import copy
from .plexception import *
from .integrate import conintegrate
from .descend import *
from .contour import *


class plmodel:
    """Perform Picard-Lefschetz integration for a given oscillatory integrand and initial contour in C^ndim.

    A plmodel object takes an oscillatory integrand of the form exp(iS) and an initial contour, and has
    methods to deform the initial contour according towards the Lefschetz thimbles (contours along which
    the integral is no longer oscillatory). Plmodel also has methods to then integrate the integrand
    over the deformed contour.

    Attributes
    ----------
        contour: plind.contour.core.contour
            The current contour along which the integration is defined. See the documentation in
            plind.contour for more information.
        expfun: function
            The function in the exponential of the oscillatory integrand. E.g. if the integrand
            is of the form exp(iS), then expfun = iS. expfun should take as its first argument a
            complex vector, z, in C^ndim (that is z is an ndarray<complex> of length ndim). It may take
            any number of additional arguments: expfun(z, *expargs).
        grad: function
            The complex gradient of the Morse function, h = Real(iS). That is grad = dh/dz. The
            gradient is required to determine the direction in which to deform the contour. The
            arguments of grad should be the same as the arguments as expfun: grad(z, *expargs).
        expargs: array_like
            Additional arguments to be passed to expfun and grad
        ndim: int
            The dimension of the complex space, i.e. C^ndim.
        trajectory: array<plind.contour.core.contour>
            A list of contours. When plmodel.descend() is called, the contour at each time step is
            appended to trajectory so that the deformation history of the contour is tracked.
        integral: tuple<complex>
            Tuple containing the value of the integral and an associated error. Default is None. When plmodel.integrate() is called integral is updated
            to be the value of the integral at the current contour.
        dt: float
            Default is None, but is updated when plmodel.descend() is called to be the current timestep
            parameter for the descent algorithm.
        delta: float
            Default is None, but is updated when plmodel.descend() is called to be the delta parameter
            used in calling that function.
        thresh: float
            Default is None, but is updated when plmodel.descend() is called to be the thresh parameter
            used in calling that function.

    """

    def __init__(self, contour, expfun, grad=None, expargs=[]):
        self.contour = contour
        self.expfun = expfun
        self.grad = grad
        self.expargs = expargs
        if np.shape(self.expargs) == ():
            self.expargs = [self.expargs]
        self.ndim = contour.simplices.shape[1]-1  # ndim = number of vertices in simplex minus 1
        self.trajectory = np.array([copy(contour)])
        self.integral = None
        # Parameters used in last plmodel.descend() call
        self.dt = None
        self.delta = None
        self.thresh = None

    # Simple functions for retrieving attributes
    def get_contour(self):
        return self.contour

    def get_expfun(self):
        return self.expfun

    def get_grad(self):
        return self.grad

    def get_trajectory(self):
        return self.trajectory

    def get_integral(self):
        return self.integral

    def get_intfun(self):
        """Return integrand function, i.e. np.exp(self.expfun(z, *self.expargs)).

        Returns
        -------
            intfun: function
                intfun = np.exp(expfun)

        """
        def intfun(z, *args):
            return np.exp(self.expfun(z, *args))
        return intfun

    def get_morse(self):
        """Return the morse function, i.e. the real part of expfun.

        Returns
        -------
            morse: function
                morse = np.real(expfun)

        """
        def morse(z, *args):
            return np.real(self.expfun(z, *args))
        return morse

    def get_morsevals(self):
        h = self.get_morse()
        return h(self.contour.points.T, *self.expargs).flatten()

    # Functions for performing the PL integration
    def descend(self, delta, thresh, tmax, dt_init, verbose=True):
        """Deform the contour according to the Picard-Lefschetz rule (flow the points along the gradient of the Morse function).

        Upon calling plmodel.descend(), the contour is deformed along the gradient of the Morse function with an adaptive mesh
        algorithm to prevent points in the contour from becoming too sparse. At each time step, plmodel.contour is updated to
        be the current contour, and the current contour is appended to plmodel.trajectory.

        Parameters
        ----------
            delta: float
                The maximum Euclidean distance points joined by an edge in the contour are allowed to be from each other.
            thresh: float
                The minimum value of expfun(z, *expargs) that a point z in the contour is allowed to be.
            tmax: float
                The time to integrate the flow equation to.
            dt_init: float
                The initial time step for the deformation. For non-adaptive timestep methods, dt is dt_init at all time during the flow, and total number of steps is celing(tmax/dt_init).
            verbose: bool
                If True, will print the total steps/time taken to descend.
        """
        h = self.get_morse()  # get the Morse function, h = real(expfun)
        gradh = self.grad

        # Descend according to the Picard-Lefschetz rule for Nstep time steps
        t = 0
        i = 0
        dt = dt_init
        while t <= tmax:
            #stop at tmax
            dt = min(dt, tmax-dt)

            # perform stepping
            self.contour.points, dt = flow(self.contour.points, gradh, dt, expargs=self.expargs)  # perform pushing

            # remove points from the contour for which self.expfun evaluated at those points
            # is below the threshold
            if self.ndim == 1:
                hval = h(self.contour.points[:,None], *self.expargs)
            else:
                hval = h(self.contour.points.T, *self.expargs)
            bad_points = np.where(hval < thresh)[0]  # find the points to remove
            if len(bad_points) > 0:
                self.contour.remove_points(bad_points)  # remove points

            # perform the adaptive mesh refinement: simplices in the contour with edge lengths greater than
            # delta are split in half
            self.contour.refine_edges(delta)

            # add new contour to trajectory
            self.trajectory = np.append(self.trajectory, copy(self.contour))

            t += dt
            i += 1
        Nstep = i
        # Store descend parameters
        self.dt = dt
        self.thresh = thresh
        self.delta = delta
        if verbose:
            print('total steps:', Nstep, 'current time:', t)


    def integrate(self, intfun=None):
        """Perform contour integration along the current contour, using the Grundmann-Moeller
        method for integration of simplices.

        When called, plmodel.integrate() performs the contour integration and updates plmodel.integral to be
        the result of the integration.

        Parameters
        ----------
            intfun: function, optional
                Function of the form intfun(z, *self.expargs) to integrate over. If None is given,
                then intfun = plmodel.get_intfun(). The default is None.

        """
        if intfun is None:
            self.intfun = self.get_intfun()
        else:
            self.intfun = intfun

        integral, gm_err = conintegrate(self.intfun, self.contour, args=self.expargs)

        # Estimate error from having too large a value for thresh
        if self.thresh != None:
            self.descend(self.delta, 10*self.thresh, 2*self.dt, self.dt, verbose=False)
            integral1 = conintegrate(self.intfun, self.contour, args=self.expargs)[0]
            thresh_err = 2*np.abs(integral-integral1)

            self.integral = (integral, np.sqrt(thresh_err**2 + gm_err**2))
        else:
            self.integral = (integral, gm_err)
