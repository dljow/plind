import numpy as np
from __future__ import division

class solution:

    def __init__(self, data):
        self.data = data

# Useful things to get from the bunch object solution, plus a wrapper for the contour
    def get_timepts(self):
        """Return the time points at which the integrator evolved the gradient."""
        return self.data.t

    def get_trajectory(self):
        """Return the trajectory of the gradient descent at the time points specified by self.get_timepts()."""
        # convert trajectory returned by ivp into complex form
        trajectory = self.data.y
        npts = (trajectory.shape[0]//2)
        trajectory = np.moveaxis(trajectory, -1, 0)
        trajectory = trajectory[:, :npts] + 1j * trajectory[:, npts:2*npts]
        return trajectory

    def get_contour(self):
        """Return the last path of the trajectory (ie. hopefully the Lefshetz thimble)."""
        trajectory = self.get_trajectory()
        return trajectory[-1, :]

    def get_ode_sol(self):
        """Return the OdeSolution instance provided by the solver."""
        return self.data.sol

# Questionably useful stuff from the bunch object below this line:
    def get_t_events(self):
        """Return the t_events from the solver, that is, for each event type a list of arrays at which an event of that type was detected."""
        return self.data.t_events

    def get_nfev(self):
        """Return the number of evaluations of the right-hand side."""
        return self.data.nfev


    def get_njacobev(self):
        """Return the number of evaluations of the Jacobian."""
        return self.data.njev

    def get_nlu(self):
        """Return number of LU decompositions."""
        return self.data.nlu

    def get_status(self):
        """Return the status of the solver."""
        return self.data.status


    def sget_message(self):
        """Return the solver message."""
        return self.data.message

    def get_success(self):
        """Return a boolean describing if the solver succeeded.s"""
        return self.data.success
