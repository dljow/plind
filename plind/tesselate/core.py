#!/usr/bin/python3
# tesselate.py
# use the Delauney tesselation to get a 2D triangulation of the plane
import numpy as np
from scipy.spatial import Delaunay

eps = np.finfo(float).eps


def get_tesselation(points):
     """ Returns a Delaunay tesselation of a set of points, in the form of 
         simplices referring to the indices of the set points, and edges with the 
         same indices. The edges have no repeats."""


      # Use Delaunay tesselation to get points
      tri = Delaunay(points)
      simplices= tri.simplices

      # Extract all possible edge pairs from tesselation
      edges= np.append(np.append(simplices[:,[0,1]],simplices[:,[1,2]]),simplices[:,[2,0]])
      edges= np.reshape(edges, [int(edges.size/2), 2]) # reshape to pairs
      edges.sort(axis=1) # put all pairs in ascending order
      edges= np.unique(edges, axis=0) # remove duplicates
      return simplices, edges
