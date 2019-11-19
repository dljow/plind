import numpy as np
from scipy.spatial import Delaunay

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
        edges = np.reshape(edges, [int(edges.size/2), 2])  # reshape to pairs
        edges.sort(axis=1)  # put all pairs in ascending order
        edges = np.unique(edges, axis=0)  # remove duplicates

        self.points = points
        self.edges = edges
        self.simplices = simplices

    # Function to compute edge lengths
    def get_edgelengths(self):
        return np.norm(self.points[edges][:,0] - self.points[edges][:,1])

    # Function to split edges in half
    def split_edges(self, bad_edges):
        for edge in bad_edges:
           p0, p1= self.points[edge]
           mid=(p0+p1)/2
           simplices_tochange= []

           # for loop can probably be made more efficient
           for simplex in self.simplices:
              if (edge[0] in simplex) and (edge[1] in simplex):
                simplices_to_change.append(simplex)


    # Function to remove points
    def remove_points(self, bad_points):
        bad_edges = np.unique(np.array(np.where(np.isin(self.edges, bad_points)))[0])
        bad_simplices = np.unique(np.array(np.where(np.isin(self.simplices, bad_points)))[0])
        self.points = np.delete(self.points, bad_points)
        self.edges = np.delete(self.edges, bad_edges)
        self.simplices = np.delete(self.simplices, bad_simplices)
        pass

    # Function to refine edges
    def refine_edges(self, delta):
        # Add points to the points that are too far away
        lengths = self.get_edgelengths()
        bad_edges= edges[lengths > delta]
        self.split_edges(bad_edges)

        # identify poles, remove
        bad_points = []
        self.remove_points(bad_points)
        # (!!!) add here
