import numpy as np
from scipy.spatial import Delaunay
from math import factorial as fact
import itertools
from ..helpers.core import unordered_pairing

class contour:
    """A contour (surface) in C^ndim for the purposes of gradient flow and integration.

           Attributes
           ----------
           points: numpy.complex64 array
                  Array of points in N-dimensional complex space, where each point is a numpy array of the
                  appropriate dimension.

           edges: numpy.int array
                  (M, (ndim+1) choose 2, 2)-array, where each (ndim+1 choose 2, 2)-subarray is a list of edges,
                  where each edge is a tuple of integers referring to the indices in self.points of the points
                  comprising the edge

           simplices: numpy.int array
                  (M, (ndim+1))-array, where each element is a list of integers referring to the indices in
                  self.points of the points comprising a simplex. There are M simplices.

           ndim: int
                  The dimension of the complex space, ie. C^ndim.
    """

    def __init__(self, points=np.array([]), simplices=np.array([[]])):
        self.points = points
        self.edges = np.array([list(itertools.combinations(simp, 2)) for simp in simplices])
        self.simplices = simplices
        self.ndim = simplices.shape[1]-1

    def init_contour(self, points):
        """Initialize a contour object from a set of points only, using a Delaunay triangulation to generate edges
           and simplices.

           Parameters
           ----------

           points: numpy.complex64 array
                  Array of points in N-dimensional complex space, where each point is a numpy array of the
                  appropriate dimension.

           Returns
           -------
           plind.contour Object


           See Also
           --------

            __init__: Initializes a plind.contour given points, edges, and simplices.

        """
        # Use Delaunay tesselation to get points
        tri = Delaunay(points)
        simplices = tri.simplices

        # Assign all quantities
        self.points = points
        self.edges = np.array([list(itertools.combinations(simp, 2)) for simp in simplices])
        self.simplices = simplices
        self.ndim = np.shape(simplices)[1]-1

    def get_edgelengths(self):
        """Return the edge lengths of the contour, as a numpy array.

           Returns
           -------

           norm_diff: np.float64 array
                (M, ndim+1 choose 2)-array, where each element is a list of (ndim+1 choose 2) values
                corresponding to the lengths of the edges in each simplex

        """
        diff = self.points[self.edges[:, :, 0]] - self.points[self.edges[:, :, 1]]
        if self.ndim == 1:
            lengths = np.abs(np.sqrt(diff**2))
        else:
            lengths = np.abs(np.sqrt(np.sum([diff[:, :, i]**2 for i in np.arange(diff.shape[2])], axis=0)))
        return lengths

    # this bit of code is taken directly from quadpy. How does credit work here?
    def get_vols(self, imag=False):
        """Compute the volume of the simplices via the Cayley-Menger determinant

           Parameters
           ----------

           imag: bool
               If False, the real-valued volume will be computed for simplices embedded in R^(2*ndim).
               If True, the volume will be imaginary and will be computed for simplices embedded in C^(ndim)

        Returns
        -------

           vol: np.float array / np.complex array
                List of volumes associated with each simplex in self.simplices. If imag==True the volumes
                will be complex valued.
        """
        simps = np.stack(self.points[self.simplices], axis=-2)
        simplex = np.asarray(simps)

        # compute all edge lengths
        edges = np.subtract(simplex[:, None], simplex[None, :])
        ei_dot_ej = np.einsum("...k,...k->...", edges, edges)

        j = simplex.shape[0] - 1
        a = np.empty((j + 2, j + 2) + ei_dot_ej.shape[2:], dtype=complex)
        a[1:, 1:] = ei_dot_ej
        a[0, 1:] = 1.0
        a[1:, 0] = 1.0
        a[0, 0] = 0.0

        a = np.moveaxis(a, (0, 1), (-2, -1))
        det = np.linalg.det(a)

        vol = np.sqrt((-1.0) ** (j + 1) / 2 ** j / fact(j) ** 2 * det)
        if imag:
            return vol
        else:
            return np.abs(vol)

    def rm_reindex(self, arr, bad_point_ind):
        """Given an array of simplices or edges, re-index the array based on the removal of points. This is
           a linear shifting of indices.

           Parameters
           ----------
           arr : np.int64 array
                 Array of point indices, corresponding to either simplices or edges
           bad_point_ind : np.int64 array
                 Array of point indices to be removed

           Returns
           -------

           arr: np.int64 array
                 The original arr, re-indexed to have the points removed.

           See Also
           --------
           remove_points: Removes the points at the indices bad_point_ind from the contour.

        """
        # make a set for speed purposes
        bad_point_ind_set = np.unique(bad_point_ind)
        amount_to_subtract = []
        counter = 0

        # loop once of all the possible indices in arr, if there is a bad index,
        # i know which position it's at, and how many bad indices there are below it
        for val in range(np.max(arr)+1):
            if val in bad_point_ind_set:
                counter += 1
            amount_to_subtract.append(counter)
        amount_to_subtract = np.array(amount_to_subtract)

        # arr then indexes amount_to_subtract for how much to subtract
        arr = arr - amount_to_subtract[arr]
        return arr

    def remove_points(self, bad_point_ind):
        """Given indices of points to be removed, remove them from self.points.

           Parameters
           ----------
           bad_point_ind : np.int64 array
                 Array of point indices to be removed

           See Also
           --------
           refine_edges: Refines the contour by splitting the edges to be finer, and removing points above
                         a certain threshold.

        """
        #bad_edge_ind = np.unique(np.array(np.where(np.isin(self.edges, bad_point_ind)))[0])
        bad_simp_ind = np.unique(np.array(np.where(np.isin(self.simplices, bad_point_ind)))[0])
        
        # remove bad edges and simplices
        self.points = np.delete(self.points, bad_point_ind, axis=0)
        self.edges = np.delete(self.edges, bad_simp_ind, axis=0)
        self.simplices = np.delete(self.simplices, bad_simp_ind, axis=0)
        if (len(self.edges) == 0): 
            raise RuntimeError("All domain edges removed. Try setting your threshold to a lower value.")
        
        # Set the indices of the points accordingly
        self.edges = self.rm_reindex(self.edges, bad_point_ind)
        self.simplices = self.rm_reindex(self.simplices, bad_point_ind)

    # Function to refine edges
    def refine_edges(self, delta):
        """Refines the contour by splitting the edges to be finer, and removing points above
                         a certain threshold.

           Parameters
           ----------
           delta: np.float64
                 Threshold size over which edges are split.

        """
        # construct edges from simplices and compute lengths
        if self.ndim == 1:
            self.points = self.points.flatten()

        diff = self.points[self.edges[:, :, 0]] - self.points[self.edges[:, :, 1]]
        if self.ndim == 1:
            lengths = np.abs(np.sqrt(diff**2))
        else:
            lengths = np.abs(np.sqrt(np.sum([diff[:, :, i]**2 for i in np.arange(diff.shape[2])], axis=0)))

        # find edges that exceed delta
        where = np.where(lengths > delta)
        if len(where[0]) > 0:
            all_bad_edges = self.edges[where]
            all_edge_key = unordered_pairing(all_bad_edges[:, 0], all_bad_edges[:, 1])  # unique indetifier for each edge

            # make it so there is only one edge per simplex, and remove secondary edges as flagged edges
            bad_simp_ind, prim_edge_ind = np.unique(where[0], return_index=True)
            dupl_edge_key = np.delete(all_edge_key, prim_edge_ind)
            dupl_edge_ind = np.where(np.in1d(all_edge_key, dupl_edge_key))[0]

            bad_simp_ind = np.delete(where[0], dupl_edge_ind)
            bad_simps = self.simplices[bad_simp_ind]
            bad_edges = np.delete(all_bad_edges, dupl_edge_ind, axis=0)
            edge_key = np.delete(all_edge_key, dupl_edge_ind)

            # create new points to add to simplices and also determine the index of these new points once added to self.points
            uni_edge_key, unique_inds = np.unique(edge_key, return_index=True)  # identify unique edges (as edges may be shared by simplices)
            uni_bad_edges = bad_edges[unique_inds]

            new_pnts = (self.points[uni_bad_edges[:, 0]] + self.points[uni_bad_edges[:, 1]])/2
            newpnts_inds = np.unique(edge_key,return_inverse=True)[1] + len(self.points)

            A = np.hstack([newpnts_inds[:, None], bad_simps[np.where(bad_simps != bad_edges[:, 0][:, None])].reshape(bad_simps.shape[0], bad_simps.shape[1]-1)])
            B = np.hstack([newpnts_inds[:, None], bad_simps[np.where(bad_simps != bad_edges[:, 1][:, None])].reshape(bad_simps.shape[0], bad_simps.shape[1]-1)])

            new_simps = np.concatenate([A, B])
            new_edges = np.rollaxis(new_simps.T[np.transpose(np.triu_indices(self.ndim+1, 1))], -1) # this line is a bit of black magic, it's from https://stackoverflow.com/questions/16003217/n-d-version-of-itertools-combinations-in-numpy

            self.simplices = np.concatenate([np.delete(self.simplices, bad_simp_ind, axis=0), new_simps])
            self.points = np.concatenate([self.points, new_pnts])
            self.edges = np.concatenate([np.delete(self.edges, bad_simp_ind, axis=0), new_edges])


            if self.ndim == 1:
                self.points = self.points[:, None]

        else:

            if self.ndim == 1:
                self.points = self.points[:, None]
