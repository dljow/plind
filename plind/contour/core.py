import numpy as np
from scipy.spatial import Delaunay
import array
from numba import jit
from ..numba_utils import *


@jit(nopython=True)
def split_edges_nb(points, edges, simplices, ndim, bad_edges, indicies):
    """Numba compiled function.
    For an array of bad_edges flagged as too long for a given mesh spacing, split the edges in half
          and modify the simplices, points arrays accordingly.

          Note that this function deals with the simplices one at a time to avoid geometric conflicts, so
          if two edges within the same simplex are flagged to be split, only one of them will be split (the
          first one to appear in the list). The assumption is that the other wedge will be flagged and split
          on the subsequent time-step.

       Parameters
       ----------
       bad_edges: np.int64 array
             Array of the edges (point indices in pairs) that are too long.
       indices : np.int64 array
             Array of edge indices to be removed.


       See Also
       --------
       refine_edges: Refines the contour by splitting the edges to be finer, and removing points above
                     a certain threshold.

    """

    used_simps = np.zeros(1, np.int64) - 1
    uni_bad_edges = np.zeros((1, 2), np.int64)
    uni_bad_simps = np.zeros((1, ndim+1), np.int64)

    # The surface becomes inconsistent if two edges within the same simplex are flagged at once.
    # To combat this, only the first edge is split, and the second is assumed to be caught
    # by the subsequent time step
    for i in np.arange(bad_edges.shape[0]):
        bad_edge = bad_edges[i]
        # Keep track of the simplices associated with an edge
        simplices_tag = np.where(sum_nb(isin_nb(simplices, bad_edge)) > 1)[0]

        if not np.any(in1d_nb(simplices_tag, used_simps)):  # Check if we have used this simplex
            # Flag simplices_tag to not reuse simplices
            if i == 0:
                used_simps = simplices_tag
            else:
                used_simps = np.concatenate((used_simps, simplices_tag))
            # Add edge to new bad edge array
            if i == 0:
                uni_bad_edges[0] = bad_edge
            else:
                uni_bad_edges = np.concatenate((uni_bad_edges.ravel(), bad_edge)).reshape((uni_bad_edges.shape[0]+1, 2))
            # Remove bad edges from the list of edges
            edges_tag = sum_nb(isin_nb(edges, bad_edge)) == 2
            edges = edges[~(edges_tag)]

            # Add simplice(s) with the proper extras populated
            if i == 0:
                uni_bad_simps[0] = simplices[simplices_tag[0]]
                for j in np.arange(1, len(simplices_tag)):
                    uni_bad_simps = np.concatenate((uni_bad_simps.ravel(), simplices[simplices_tag[j]])).reshape((uni_bad_simps.shape[0]+1, uni_bad_simps.shape[1]))
            else:
                for tag in simplices_tag:
                    uni_bad_simps = np.concatenate((uni_bad_simps.ravel(), simplices[tag])).reshape((uni_bad_simps.shape[0]+1, uni_bad_simps.shape[1]))

    # The only edges this function treats and removes are uni_bad_edges,
    # that is, the edges that are unique with respect to their simplices.
    uni_bad_simps = uni_bad_simps.reshape(-1, ndim+1)
    uni_bad_edges = uni_bad_edges.reshape(-1, 2)

    # Add the midpoints of all the bad edges
    midpts_ind = np.arange(points.shape[0], points.shape[0] + uni_bad_edges.shape[0], 1)
    midpts = (points[uni_bad_edges[:, 0]] + points[uni_bad_edges[:, 1]])/2
    points = np.concatenate((points.ravel(), midpts.ravel())).reshape(len(points)+len(midpts), ndim)

    # Add the two edges that replace all the bad edges
    edges_1 = np.concatenate((midpts_ind, uni_bad_edges[:, 0])).reshape(2, -1).T
    edges_2 = np.concatenate((midpts_ind, uni_bad_edges[:, 1])).reshape(2, -1).T
    edges = np.concatenate((edges, edges_1, edges_2), axis=0)

    # Delete all the bad edges
    simplices = delete_row_nb(simplices, used_simps)


    # vertices which are not part of the bad edges in the bad simplices
    for i in np.arange(0, len(uni_bad_edges)):
        # Get all points in the simplex not associated to the edge ("outliers")
        bad_edge = uni_bad_edges[i]

        outliers = bool_index_nb(uni_bad_simps, isin_nb(uni_bad_simps, bad_edge, invert=True) *
                                 rowarr_transpose_nb(sum_nb(isin_nb(uni_bad_simps, bad_edge, invert=True)) == ndim+1-2))
        uni_outliers = np.unique(outliers)

        # ndim - 1 outliers will exist in every edge
        if ndim == 1:
            num_simps = 1
        else:
            num_simps = int(outliers.size/(ndim-1))
        num_outliers = uni_outliers.size

        outliers_T = outliers.reshape(num_simps, ndim-1)
        uni_outliers_T = uni_outliers.reshape(num_outliers, 1)

        # add new edge for every outlier
        edges_outliers = np.concatenate((uni_outliers_T, midpts_ind[i]*np.ones((num_outliers, 1), dtype=np.int64)), axis=1)
        edges = np.concatenate((edges, edges_outliers), axis=0)

        # Simplices per outlier row
        simp_1 = np.concatenate((midpts_ind[i]*np.ones((num_simps, 1), dtype=np.int64), outliers_T, bad_edge[0]*np.ones((num_simps, 1), dtype=np.int64)), axis=1)
        simp_2 = np.concatenate((midpts_ind[i]*np.ones((num_simps, 1), dtype=np.int64), outliers_T, bad_edge[1]*np.ones((num_simps, 1), dtype=np.int64)), axis=1)

        simplices = np.concatenate((simplices, simp_1, simp_2), axis=0)
    return points, edges, simplices


class contour:
    """A contour (surface) in C^ndim for the purposes of gradient flow and integration.

           Attributes
           ----------
           points: numpy.complex64 array
                  Array of points in N-dimensional complex space, where each point is a numpy array of the
                  appropriate dimension.

           edges: numpy.int16 array
                  Array of two-element arrays, with integer elements that refer to the points contained in
                  an edge.

           simplices: numpy.int16 array
                  Array of ndim-element arrays, with integer elements that refer to the points contained in
                  an edge.
           ndim: int
                  The dimension of the complex space, ie. C^ndim.
    """

    def __init__(self, points=np.array([]), edges=np.array([[]]), simplices=np.array([[]])):
        self.points = points
        self.edges = edges
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

        # Extract all possible edge pairs from tesselation
        edges = np.reshape(np.append(np.append(simplices[:, [0, 1]], simplices[:, [1, 2]]), simplices[:, [2, 0]]), [int(edges.size/2), 2])
        edges.sort(axis=1)  # put all pairs in ascending order
        edges = np.unique(edges, axis=0)  # remove duplicates

        # Assign all quantities
        self.points = points
        self.edges = edges
        self.simplices = simplices
        self.ndim = np.shape(simplices)[1]-1

    def get_edgelengths(self):
        """Return the edge lengths of the contour, as a numpy array.

           Returns
           -------

           norm_diff: np.float64 array
                     Numpy array of the lengths of the edges.

           See Also
           --------
           refine_edges: Refines the edges of a plind.contour to be smaller than a given delta.

        """
        differences = (self.points[self.edges][:, 0] - self.points[self.edges][:, 1])

        if self.ndim == 1:
            norm_diff= np.sqrt(differences**2)
        else:
            norm_diff= np.sqrt(np.sum([differences[:, i]**2 for i in np.arange(0, self.ndim)], 0))
        return norm_diff

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
        arr = arr - np.count_nonzero(arr.ravel()[:, np.newaxis] > bad_point_ind,axis=1).reshape(arr.shape)
        return arr

    def split_edges(self, bad_edges, indices):
        new_points, new_edges, new_simplices = split_edges_nb(self.points, self.edges, self.simplices, self.ndim, bad_edges, indices)
        self.points = new_points
        self.edges = new_edges
        self.simplices = new_simplices

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
        bad_edge_ind = np.unique(np.array(np.where(np.isin(self.edges, bad_point_ind)))[0])
        bad_simp_ind = np.unique(np.array(np.where(np.isin(self.simplices, bad_point_ind)))[0])
        # remove bad edges and simplices
        self.points = np.delete(self.points, bad_point_ind, axis=0)
        self.edges = np.delete(self.edges, bad_edge_ind, axis=0)
        self.simplices = np.delete(self.simplices, bad_simp_ind, axis=0)
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
        # Add points to the points that are too far away
        lengths = self.get_edgelengths()
        bad_edges = self.edges[lengths > delta]
        if len(bad_edges) > 0:
            self.split_edges(bad_edges, np.where(lengths > delta))

        # identify poles, remove
        bad_points = []
        self.remove_points(bad_points)
