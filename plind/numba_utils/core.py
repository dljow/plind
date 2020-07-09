import numpy as np
from numba import jit

@jit(nopython=True)
def isin_nb(simplices, edge, invert=False):
    print('hi')
    isin_arr = [False]
    for a in simplices.flatten():
        bool_val = False
        for e in edge:
            if a == e:
                bool_val = True
        isin_arr.append(bool_val)

    isin_arr = np.array(isin_arr[1:])
    isin_arr = isin_arr.reshape(simplices.shape)
    if invert:
        return np.invert(isin_arr)
    else:
        return isin_arr

# def isin_nb(simplices, edge, invert=False):
#     isin_arr = [False]
#     for i in np.arange(simplices.shape[0]):
#         simp = simplices[i]
#         for a in simp:
#             bool_val = False
#             for e in edge:
#                 if a == e:
#                     bool_val = True
#             isin_arr.append(bool_val)
#
#     isin_arr = np.array(isin_arr[1:])
#     isin_arr = isin_arr.reshape(simplices.shape)
#     if invert:
#         return np.invert(isin_arr)
#     else:
#         return isin_arr

@jit(nopython=True)
def sum_nb(isin_arr):
    summed = np.zeros(isin_arr.shape[0], dtype=np.int64)
    for i in np.arange(0, isin_arr.shape[0]):
        summed[i] = np.sum(isin_arr[i])
    return summed

@jit(nopython=True)
def in1d_nb(arr0, arr1):
    out_arr = np.zeros(len(arr0), np.bool_)
    for i in np.arange(0, len(arr0)):
        for j in np.arange(0, len(arr1)):
            if arr0[i] == arr1[j]:
                out_arr[i] = True
    return out_arr

@jit(nopython=True)
def delete_row_nb(arr, num):
    mask = np.zeros(arr.shape[0], dtype=np.int64) == 0
    mask[num] = False
    return arr[mask]

@jit(nopython=True)
def rowarr_transpose_nb(rowarr):
    colarr = np.zeros((len(rowarr), 1), rowarr.dtype)
    for i in np.arange(0, len(rowarr)):
        colarr[i, 0] = rowarr[i]
    return colarr

@jit(nopython=True)
def bool_index_nb(arr, bool_arr):
    out = np.zeros(np.sum(bool_arr.flatten()), arr.dtype)
    j = 0
    for i in np.arange(0, bool_arr.size):
        if bool_arr.flatten()[i]:
            out[j] = arr.flatten()[i]
            j += 1
    return out
