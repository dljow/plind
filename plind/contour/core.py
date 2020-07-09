import numpy as np
from scipy.spatial import Delaunay
import array
from .split_edges import split_edges_nb as split_edges_nb


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
        new_points, new_edges, new_simplices = split_edges_nb(self.points, self.edges.astype(np.int64),
                                                           self.simplices.astype(np.int64), self.ndim,
                                                           bad_edges.astype(np.int64), indices.astype(np.int64))
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
            self.split_edges(bad_edges, np.where(lengths > delta)[0])

        # identify poles, remove
        bad_points = []
        self.remove_points(bad_points)
