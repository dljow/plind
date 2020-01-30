import numpy as np
from scipy.spatial import Delaunay
from time import time

class contour:

    def __init__(self, points=np.array([]), edges=np.array([[]]), simplices=np.array([[]])):
        self.points = points
        self.edges = edges
        self.simplices = simplices
        self.ndim = 0 

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
        self.ndim = np.shape(simplices)[1]-1

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
        used_simps = np.array([],dtype=np.int)
        uni_bad_edges=np.array([],dtype=np.int)
        uni_bad_simps=np.array([],dtype=np.int)
        for i, bad_edge in enumerate(bad_edges):
            simplices_tag = np.count_nonzero(np.isin(self.simplices, bad_edge), axis=-1) > 1
            simplices_tag = np.where(simplices_tag)[0]
            if not np.any(np.in1d(simplices_tag,used_simps)): # Check if we have used this simplex
                # Flag simplices_tag to not reuse simplices
                used_simps = np.append(used_simps, simplices_tag)
                # Add edge to new bad edge array
                uni_bad_edges = np.append(uni_bad_edges, bad_edge)
                # Also, we want to delete bad_edges from the original list of edges,
                edges_tag = np.count_nonzero(np.isin(self.edges, bad_edge),axis=-1) == 2
                edges_tag = np.where(edges_tag)[0]
                self.edges = np.delete(self.edges, edges_tag, axis=0)
                # Add simplice(s) with the proper extras populated
                for j in range(self.ndim):
                    if np.size(simplices_tag)>j:
                        uni_bad_simps = np.append(uni_bad_simps,self.simplices[simplices_tag[j]], axis=0)
                    else:
                        junk_simp=np.full(np.shape(self.simplices[0]),np.nan)
                        junk_simp[0]=bad_edge[0]
                        junk_simp[1]=bad_edge[1]
                        uni_bad_simps = np.append(uni_bad_simps, junk_simp, axis=0)
        uni_bad_simps = uni_bad_simps.reshape(-1,3)
        uni_bad_edges = uni_bad_edges.reshape(-1,2)

        # add points 
        midpts_ind= np.arange(np.shape(self.points)[0], np.shape(self.points)[0]+np.shape(uni_bad_edges)[0],1)
        midpts = (self.points[bad_edges[:,0]] + self.points[bad_edges[:,1]])/2
        self.points=np.append(self.points, midpts, axis=0)

        # add edges
        edges_1 = np.sort(np.append(midpts_ind, uni_bad_edges[:,0],axis=0),axis=0)
        edges_1 = np.reshape(edges_1, [np.size(midpts_ind),2])

        edges_2 = np.sort(np.append(midpts_ind, uni_bad_edges[:,1],axis=0),axis=0)
        edges_2= np.reshape(edges_2, [np.size(midpts_ind), 2])
        self.edges= np.append(self.edges, edges_1,axis=0)
        self.edges= np.append(self.edges, edges_2,axis=0)

        # used_simps conveniently tracks bad simplices after unifiquation
        self.simplices = np.delete(self.simplices, used_simps, axis=0)

        # vertices which are not part of the bad edges in the bad simplices
        for i, bad_edge in enumerate(uni_bad_edges):
                # get all outliers for all simplices associated to the bad edge
                outliers= uni_bad_simps[np.isin(uni_bad_simps, bad_edge, invert=True) *
                        (np.count_nonzero(np.isin(uni_bad_simps, bad_edge, invert=True),axis=-1)==self.ndim+1-2)[:,np.newaxis]]
                # ndim- 1 outliers will exist in every edge
                num_simps = int(np.size(outliers)/(self.ndim-1))
                
                # add new edge for every outlier
                edges_outliers = np.sort(np.array(np.meshgrid(outliers, midpts_ind)).T.reshape(-1,2),axis=0)
                self.edges= np.append(self.edges, edges_outliers,axis=0)


                outliers = np.reshape(outliers, [num_simps, self.ndim-1])

                # Simplices per outlier row
                simp_1= np.sort(np.append(np.append(outliers, midpts_ind[i]*np.ones([num_simps,1]),axis=1), bad_edge[0]*np.ones([num_simps,1]),axis=1),axis=1)
                simp_2= np.sort(np.append(np.append(outliers, midpts_ind[i]*np.ones([num_simps,1]),axis=1), bad_edge[1]*np.ones([num_simps,1]),axis=1),axis=1)
               
                self.simplices=np.append(self.simplices, simp_1,axis=0)
                self.simplices= np.append(self.simplices, simp_2, axis=0)

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
            self.split_edges(bad_edges, np.where(lengths>delta))

        # identify poles, remove
        bad_points = []
        self.remove_points(bad_points)
