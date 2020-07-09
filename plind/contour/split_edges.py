import numpy as np
from numba import jit
from ..numba_utils import *
from time import time

def split_edges(points, edges, simplices, ndim, bad_edges, indicies):
    used_simps = np.array([], dtype=np.int)
    uni_bad_edges = np.array([], dtype=np.int)
    uni_bad_simps = np.empty((0, ndim+1), int)

    # The surface becomes inconsistent if two edges within the same simplex are flagged at once.
    # To combat this, only the first edge is split, and the second is assumed to be caught
    # by the subsequent time step
    for i, bad_edge in enumerate(bad_edges):
        # Keep track of the simplices associated with an edge
        simplices_tag = np.isin(simplices, bad_edge).sum(axis=-1) > 1
        simplices_tag = np.where(simplices_tag)[0]

        if not np.any(np.in1d(simplices_tag, used_simps)):  # Check if we have used this simplex
            # Flag simplices_tag to not reuse simplices
            used_simps = np.append(used_simps, simplices_tag)
            # Add edge to new bad edge array
            uni_bad_edges = np.append(uni_bad_edges, bad_edge)
            # Remove bad edges from the list of edges
            edges_tag = np.isin(edges, bad_edge).sum(axis=-1) == 2
            edges = edges[~(edges_tag)]

            # Add simplice(s) with the proper extras populated
            uni_bad_simps = np.append(uni_bad_simps, simplices[simplices_tag], axis=0)

    # The only edges this function treats and removes are uni_bad_edges,
    # that is, the edges that are unique with respect to their simplices.
    uni_bad_simps = uni_bad_simps.reshape(-1, ndim+1)
    uni_bad_edges = uni_bad_edges.reshape(-1, 2)

    # Add the midpoints of all the bad edges
    midpts_ind = np.arange(np.shape(points)[0], np.shape(points)[0]+np.shape(uni_bad_edges)[0], 1, dtype=np.int)
    midpts = (points[uni_bad_edges[:, 0]] + points[uni_bad_edges[:, 1]])/2
    points = np.append(points, midpts, axis=0)

    # Add the two edges that replace all the bad edges
    edges_1 = np.sort(np.append(midpts_ind, uni_bad_edges[:, 0], axis=0).reshape(2,-1).T, axis=1)
    edges_2 = np.sort(np.append(midpts_ind, uni_bad_edges[:, 1], axis=0).reshape(2,-1).T, axis=1)
    edges = np.concatenate((edges, edges_1, edges_2), axis=0)

    # Delete all the bad edges
    simplices = np.delete(simplices, used_simps, axis=0)

    # vertices which are not part of the bad edges in the bad simplices
    for i, bad_edge in enumerate(uni_bad_edges):
            # Get all points in the simplex not associated to the edge ("outliers")
            outliers = uni_bad_simps[np.isin(uni_bad_simps, bad_edge, invert=True) *
                        (np.isin(uni_bad_simps, bad_edge, invert=True).sum(axis=-1)==ndim+1-2)[:, np.newaxis]]
            uni_outliers = np.unique(outliers)

            # ndim - 1 outliers will exist in every edge
            if ndim == 1:
                num_simps = 1
            else:
                num_simps = int(np.size(outliers)/(ndim-1))
            num_outliers = np.size(uni_outliers)

            outliers = np.reshape(outliers, [num_simps, ndim-1])
            uni_outliers = np.reshape(uni_outliers, [num_outliers, 1])

            # add new edge for every outlier
            edges_outliers = np.sort(np.append(uni_outliers, midpts_ind[i]*np.ones([num_outliers, 1], dtype=np.int),axis=1), axis=0)
            edges_outliers = edges_outliers.astype(np.int)
            edges = np.append(edges, edges_outliers, axis=0)

            # Simplices per outlier row
            simp_1 = np.sort(np.concatenate(( midpts_ind[i]*np.ones([num_simps, 1], dtype=np.int),outliers,bad_edge[0]*np.ones([num_simps, 1], dtype=np.int)),axis=1), axis=1)
            simp_2 = np.sort(np.concatenate(( midpts_ind[i]*np.ones([num_simps, 1], dtype=np.int),outliers,bad_edge[1]*np.ones([num_simps, 1], dtype=np.int)),axis=1), axis=1)

            simplices = np.concatenate((simplices, simp_1, simp_2), axis=0)
    return points, edges, simplices


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
        simplices_tag = np.where(np.sum(isin_nb(simplices, bad_edge), axis=-1) > 1)[0]


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
            edges_tag = np.sum(isin_nb(edges, bad_edge), axis=-1) == 2
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
                                 rowarr_transpose_nb(np.sum(isin_nb(uni_bad_simps, bad_edge, invert=True), axis=-1) == ndim+1-2))
        uni_outliers = np.unique(outliers)

        # ndim - 1 outliers will exist in every edge
        if ndim == 1:
            num_simps = 1
        else:
            num_simps = np.int(outliers.size/(ndim-1))
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
