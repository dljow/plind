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
        for pointind in simplex:
            if pointind not in edge:
                 return pointind

    # Reindexes simplices or edges given a list of bad_points that will be removed
    def rm_reindex(self, arr, bad_points):
        arr = arr - np.array([len(np.where(i > bad_points)[0]) for i in arr.flatten()]).reshape(arr.shape)
        return arr


    # Function to split edges in half
    def split_edges(self, bad_edges, indices):
        # Store the current number of points for indexing purposes
        num_points=np.shape(self.points)[0]
        print("bad edges")
        print(bad_edges)
        print("edges before")
        print(self.edges)
        
        # delete the bad edges from the overall edge set
        self.edges= np.delete(self.edges, indices, axis=0)
        print("edges after deletion")
        print(self.edges)
        
        # add the new midpoints to the overall points set
        bad_points=self.points[bad_edges]
        bad_points=bad_points[0]
        midpoints = (bad_points[0]+bad_points[1])/2  
        self.points=np.append(self.points, [midpoints], axis=0)
        
        # locate the simplices that need changing
        bad_edges=bad_edges.tolist()
        simplices_to_change=[(simplex,edge, bad_edges.index(edge)) for simplex in self.simplices for edge in bad_edges if ((edge[0] in simplex) and (edge[1] in simplex))]
        
        # Delete bad simplices
        bad_simplices = [simplex_tuple[0] for simplex_tuple in simplices_to_change]
        self.simplices = np.delete(self.simplices, bad_simplices, axis=0)
        
        # Construct new simplices
        new_simps1= [sorted([simp_tuple[2]+num_points, self.find_peak(simp_tuple[0], simp_tuple[1]), simp_tuple[1][0]]) for simp_tuple in simplices_to_change]
        new_simps2= [sorted([simp_tuple[2]+num_points, self.find_peak(simp_tuple[0], simp_tuple[1]), simp_tuple[1][1]]) for simp_tuple in simplices_to_change]
        self.simplices=np.append(self.simplices, new_simps1, axis=0)
        self.simplices=np.append(self.simplices, new_simps2, axis=0)
        
        # Construct new edges
        new_edges1= [sorted([self.find_peak(simp_tuple[0], simp_tuple[1]), simp_tuple[2]+num_points]) for simp_tuple in simplices_to_change]
        new_edges2=[sorted([edge[0], bad_edges.index(edge)+num_points]) for edge in bad_edges]
        new_edges3=[sorted([edge[1], bad_edges.index(edge)+num_points]) for edge in bad_edges]
        
        self.edges=np.append(self.edges, new_edges1, axis=0)
        self.edges=np.append(self.edges, new_edges2, axis=0)
        self.edges=np.append(self.edges, new_edges3, axis=0)
        print("edges at end")
        print(self.edges)



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
        if len(bad_edges)>0:
          self.split_edges(bad_edges, lengths>delta)

        # identify poles, remove
        bad_points = []
        self.remove_points(bad_points)
