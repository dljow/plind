import numpy as np
from scipy.spatial import Delaunay

class contour:

    def __init__(self, points=np.array([]), edges=np.array([[]]), simplices=np.array([[]])):
        self.points = points
        self.edges = edges
        self.simplices = simplices

    # Function that initializes contour based on list of points
    def init_contour(self, points):
<<<<<<< HEAD
        """Return a Delaunay tesselation of a set of points, in the form of simplices referring to the indices of the set points, and edges with the same indices. The edges have no repeats."""
=======
        """Return a Delaunay tesselation of a set of points, in the form ofsimplices referring to the indices of the set points, and edges with the same indices. The edges have no repeats."""
>>>>>>> e51737a1a4f68e0439ab0dbd3f5f485ccc614d42
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
<<<<<<< HEAD
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
   
           

=======
        pass

    # Function to split edges in half
    def split_edges(self, bad_edges):
>>>>>>> e51737a1a4f68e0439ab0dbd3f5f485ccc614d42
        pass

    # Function to remove points
    def remove_points(self, bad_points):
        pass

    # Function to refine edges
    def refine_edges(self, delta):
<<<<<<< HEAD

        # Add points to the points that are too far away
        lengths = self.get_edgelengths()
        bad_edges= edges[lengths > delta]
        self.split_edges(bad_edges)

        # identify poles, remove
        bad_points = []
        self.remove_points(bad_points)
        # (!!!) add here
=======
        pass
>>>>>>> e51737a1a4f68e0439ab0dbd3f5f485ccc614d42
