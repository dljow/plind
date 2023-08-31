
import numpy as np
from scipy.integrate import fixed_quad
from copy import copy
from .plexception import *
from .integrate import conintegrate
from .visualize import *
from .descend import *
from .contour import *



class PLModel:
    """Perform Picard-Lefschetz integration for a given oscillatory integrand and initial contour in C^ndim.

    A PLModel object takes an oscillatory integrand of the form exp(iS) and an initial contour, and has
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
            A list of contours. When PLModel.descend() is called, the contour at each time step is
            appended to trajectory so that the deformation history of the contour is tracked.

        integral: Tuple, (float, float)
            integral(0) The value of the integral. Default is None. When PLModel.integrate(), integral is updated
            to be the value of the integral at the current contour.

            integral(1) An error estimate for the integral. See PLModel.integrate()

        dt: float
            The time step of the last run of the contour deformation.
            Set and stored after a PLModel.descend() call, initialized as None.

        delta: float
            The maximum Euclidean distance points joined by an edge in the contour are allowed to be from each other.
            Set and stored after a PLModel.descend() call, initialized as None.

        thresh: float
            The exponential threshold for the points counted for integration. A point x is discarded if (h(x) < thresh), or
            equivalently if e^(h(x)) < e^(thresh). Set and stored after a PLModel.descend() call, initialized as None.

    """

    def __init__(self, contour, expfun, grad, expargs=[]):
        self.contour = contour
        self.expfun = expfun
        self.grad = grad
        self.expargs = expargs
        if np.shape(self.expargs) == ():
            self.expargs = [self.expargs]
        self.ndim = contour.simplices.shape[1]-1  # ndim = number of vertices in simplex minus 1
        self.trajectory = np.array([copy(contour)])
        self.integral = None

        # Parameters used in last PLModel.descend() call

        self.dt = None
        self.delta = None
        self.thresh = None


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
        """Return the values of the morse function along the contour.

        Returns
        -------
            morse: function
                morse = np.real(expfun)

        """
        h = self.get_morse()
        morsevals = h(self.contour.points.T, *self.expargs).flatten()
        return morsevals

    def descend(self, delta, thresh, tmax, dt_init, verbose=True):
        """Deform the contour according to the Picard-Lefschetz rule (flow the points along the gradient of the Morse function).

        Upon calling PLModel.descend(), the contour is deformed along the gradient of the Morse function with an adaptive mesh
        algorithm to prevent points in the contour from becoming too sparse. At each time step, PLModel.contour is updated to
        be the current contour, and the current contour is appended to PLModel.trajectory.

        Parameters
        ----------
            delta: float
                The maximum Euclidean distance points joined by an edge in the contour are allowed to be from each other. Must                   be set for every run of descend(). Stored as an attribute of the PLModel object.

            thresh: float
                The exponential threshold for the points counted for integration. A point x is discarded if (h(x) < thresh), or
                equivalently if e^(h(x)) < e^(thresh). Must be set for every run of descend(). Stored as an attribute of the                     PLModel object.

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
            new_points, dt = flow(self.contour.points.T, gradh, dt, expargs=self.expargs)  # perform pushing
            self.contour.points = new_points.T
            # remove points from the contour for which self.h evaluated at those points
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

        # Store descend parameters from the last run of descend()
        self.dt = dt
        self.thresh = thresh
        self.delta = delta
        if verbose:
            print('total steps:', Nstep, 'current time:', t)


    def integrate(self, intfun=None):
        """Perform contour integration along the current contour. Defaults to a Grundmann-Moeller
        Integration scheme, with errors reported from both the integration scheme and the threshold
        for discarding points.

        To retrieve the result, call PLModel.integral.

        Parameters
        ----------
            intfun: function, optional
                Function of the form intfun(z, *self.expargs) to integrate over. If None is given,
                then intfun = PLModel.get_intfun(). The default is None.

        """
        if intfun is None:
            self.intfun = self.get_intfun()
        else:
            self.intfun = intfun

        integral, gm_err = conintegrate(self.intfun, self.contour, args=self.expargs)

        # Estimate error from having too large a value for thresh
        if self.thresh is not None:
            if self.thresh >= -20:
                self.descend(self.delta, 10*self.thresh, 2*self.dt, self.dt, verbose=False)
                integral1 = conintegrate(self.intfun, self.contour, args=self.expargs)[0]
                thresh_err = 2*np.abs(integral-integral1)

                # Add the errors from the threshold and the error from the degree of integration together
                self.integral = (integral, np.sqrt(thresh_err**2 + gm_err**2))
            else:
                # were we doing something here?
                self.integral = (integral, gm_err)
                # pass
        else:
            self.integral = (integral, gm_err)

    def visualize(self, *args, **kwargs):
        """Plots the contour at the specified step in the evolution.

        PLModel.visualize() generates a plot of the final contour. Right now, it
        only supports a 1D visualization.

        Parameters
        ----------
            step: int, optional
                Which step in the evolution to plot. Default is -1 corresponding to the final step.

            with_background: bool, optional
                Whether to plot the morse function as the background. Default is True.

            with_contour: bool, optional
                Whether to plot contours of constant Imag(iS) as a visual guide. Default is True.

            with_thresh: bool, optional
                Whether to plot the integration threshold thresh. Default is True.

        """
        if self.ndim == 1:
            plt.figure(figsize=(3.5,3.3))
            plot_1d(self, *args, **kwargs)
            plt.show()
        elif self.ndim == 2:
            print('2D viz WIP')
            pass
        else:
            print('Only 1D and 2D visualizations supported')
