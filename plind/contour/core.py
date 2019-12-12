import numpy as np
from scipy.spatial import Delaunay
from time import time

class contour:

    def __init__(self, points=np.array([]), edges=np.array([[]]), simplices=np.array([[]])):
        self.points = points
        self.edges = edges
        self.simplices = simplices

    # Function that initializes contour based on list of points
    def init_contour(self, points):
        """Return a Delaunay tesselation of a set of points, in the form of simplices referring to the indices of the set points, and edges with the same indices. The edges have no repeats"""
        # Use Delaunay tesselation to get points
        tri = Delaunay(points)
        simplices = tri.simplices

        # Extract all possible edge pairs from tesselation
        edges = np.append(np.append(simplices[:, [0, 1]], simplices[:, [1, 2]]), simplices[:, [2, 0]])
        edges = np.reshape(edges, [int(edges.size/2), 2])  # reshape to pair python list remove list of lists
        edges.sort(axis=1)  # put all pairs in ascending order
        edges = np.unique(edges, axis=0)  # remove duplicates

        self.points = points
        self.edges = edges
        self.simplices = simplices

    # Function to compute edge lengths
    def get_edgelengths(self):
        ndim = self.points.shape[1]
        differences = (self.points[self.edges][:, 0] - self.points[self.edges][:, 1])
        return np.sqrt(np.sum([differences[:, i]**2 for i in np.arange(0, ndim)], 0))


    def find_peak(self, simplex, edge):
        for point in simplex:
            if point not in edge:
                 return point

    # Reindexes simplices or edges given a list of bad_points that will be removed
    def rm_reindex(self, arr, bad_points):
        arr = arr - np.array([sum(i > k for k in bad_points) for i in arr.flatten()]).reshape(arr.shape)
        return arr


    # Function to split edges in half
    def split_edges(self, bad_edges):
        for edge in bad_edges:
            index = np.where(np.all(edge==self.edges, axis=1))
            self.edges=np.delete(self.edges, index, axis=0)

            # calculate midpoint
            p0, p1 = self.points[edge]
            mid = (p0+p1)/2
            simplices_to_change = []

            # relevant simplices:
            # (!!!) for loop can probably be made more efficient
            for simplex in self.simplices:
                if (edge[0] in simplex) and (edge[1] in simplex):
                    simplices_to_change.append(simplex)
            if (len(simplices_to_change) >2):
               print("too many simplices")
               return

            # add points
            self.points= np.append(self.points,[mid], axis=0)
            mid_ind = np.shape(self.points)[0]-1

            # add edges, simplices
            for simplex in simplices_to_change:
               # remove old simplex
               simplex.tolist()
               index= np.where(np.all(simplex==self.simplices,axis=1))
               self.simplices = np.delete(self.simplices, index, axis=0)

               # add simplices
               simp_peak= self.find_peak(simplex, edge)
               others= simplex.copy().tolist()
               others.remove(simp_peak)
               s1=[simp_peak, mid_ind, others[0]]
               s1.sort()
               s2= [simp_peak, mid_ind, others[1]]
               s2.sort()

               self.simplices=np.append(self.simplices, [s1], axis=0)
               self.simplices=np.append(self.simplices,	[s2], axis=0)

               # add edges
               self.edges= np.append(self.edges, [[simp_peak,mid_ind]], axis=0)
               self.edges= np.append(self.edges, [[edge[0], mid_ind]], axis=0)
               self.edges= np.append(self.edges, [[edge[1], mid_ind]], axis=0)



    # Function to remove points
    def remove_points(self, bad_points):
        bad_edges = np.unique(np.array(np.where(np.isin(self.edges, bad_points)))[0])
        bad_simplices = np.unique(np.array(np.where(np.isin(self.simplices, bad_points)))[0])
        # remove bad edges and simplices
        self.points = np.delete(self.points, bad_points, axis=0)
        self.edges = np.delete(self.edges, bad_edges, axis=0)
        self.simplices = np.delete(self.simplices, bad_simplices, axis=0)
        # relable points
        t0 = time()
        self.edges = self.rm_reindex(self.edges, bad_points)
        self.simplices = self.rm_reindex(self.simplices, bad_points)
        t_rind = time()-t0
        return t_rind

    # Function to refine edges
    def refine_edges(self, delta):
        # Add points to the points that are too far away
        lengths = self.get_edgelengths()
        bad_edges = self.edges[lengths > delta]
        self.split_edges(bad_edges)

        # identify poles, remove
        bad_points = []
        self.remove_points(bad_points)
