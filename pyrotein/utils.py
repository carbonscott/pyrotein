#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from operator import itemgetter
from itertools import groupby
from .atom import constant_atomlabel, constant_aminoacid_code


def bin_image(img_orig, binning = 4, mode = 1, nan_replace = 0):
    ''' Bin an image for faster display.
    '''
    Y, X = img_orig.shape

    if mode == 0:
        img_bin = []
        for i in range(0, Y, binning):
            for j in range(0, X, binning):
                sub_img = img_orig[i : min(i + binning, Y), j : min(j + binning, X)]
                if np.all(np.isnan(sub_img)):
                    img_bin(append( (i, j, nan_replace) ))
                else:
                    img_bin.append( (i, j, np.nanmean(sub_img)) )

    if mode == 1:
        img_bin = []
        for i in range(0, Y, binning):
            img_bin_y = []
            for j in range(0, X, binning):
                sub_img = img_orig[i : min(i + binning, Y), j : min(j + binning, X)]
                if np.all(np.isnan(sub_img)): 
                    img_bin_y.append( nan_replace )
                else:
                    img_bin_y.append( np.nanmean(sub_img) )
            img_bin.append(img_bin_y)

    return np.array(img_bin)


def read_file(file, numerical = False):
    '''Return all lines in the user supplied parameter file without comments.
    '''
    lines = []
    with open(file,'r') as fh:
        for line in fh.readlines():
            # Separate entries by spaces and remove commented lines...
            words = line.replace('#', ' # ').split()

            # Omit any thing coming after the pound sign in a line...
            if "#" in words: words = words[  : words.index("#")]

            # Save non-empty line...
            if numerical: words = [ float(word) for word in words ]
            if len(words) > 0: lines.append(words)

    return lines




# [[[ Matrix operation ]]]

def mat2tril(mat, keepdims = False, offset = 0):
    ''' Convert a matrix into a lower triangular matrix.
        mode:
        - False: return one-dimensional array.
        - True : return trigular matrix.
    '''
    # Convert full matrix to lower triangular matrix...
    res = mat * np.tri(len(mat), len(mat), offset)

    # Convert the lower triangular matrix into a one-dimensional array...
    if not keepdims: res = res[np.tril_indices(len(mat), offset)]

    return res




def array2tril(ary, length, offset = 0):
    ''' Convert a one-dimensional array into a lower triangular matrix.
    '''
    # Create an empty matrix with edge size of len...
    res = np.zeros((length, length))

    # Collect indices for members in lower triangular matrix with offset...
    ver_i, hor_i = np.tril_indices(length, offset)

    # Find the smaller length for valid assignment...
    capacity        = len(ver_i)
    area            = len(ary)
    rightmost_index = np.min([area, capacity])

    # Update empty matrix with values in the input array...
    res[ver_i[:rightmost_index], hor_i[:rightmost_index]] = ary[:rightmost_index]

    return res




def fill_nan_with_mean(mat, axis = 0):
    ''' Fill np.nan with mean value along `axis`.
        Support two-dimensional matrix only.
    '''
    # Assert mat is 2d...
    assert len(mat.shape) == 2, "fill_nan_with_mean ONLY supports 2D matrix."

    # Assert axis is either 0 or 1 only...
    assert axis == 0 or axis == 1, "fill_nan_with_mean ONLY allows 0 or 1 for axis."

    # Obtain the axis mean...
    axis_mean = np.nanmean(mat, axis = axis)

    # Find the indices that has values of np.nan...
    nan_i = np.where(np.isnan(mat))

    # Replace np.nan with mean...
    rep_axis = 1 - axis
    mat[nan_i] = np.take(axis_mean, nan_i[rep_axis])

    return None




def fill_nan_with_zero(mat):
    ''' Fill np.nan with zero along `axis`.
    '''
    # Assert mat is 2d...
    assert len(mat.shape) == 2, "fill_nan_with_mean ONLY supports 2D matrix."

    # Find the indices that has values of np.nan...
    nan_i = np.where(np.isnan(mat))

    # Replace np.nan with mean...
    mat[nan_i] = 0.0

    return None




def group_consecutive_integer(data):
    ''' As indicated by the function name.  Refer to 

        https://docs.python.org/2.6/library/itertools.html#examples

        for the method.  
    '''
    data_export = []
    for k, g in groupby(enumerate(data), lambda x: x[0]-x[1]):
        data_export.append( list(map(itemgetter(1), g)) )

    return data_export




def get_key_by_max_value(obj_dict):
    ''' A utility to fetch key corresponding to the max value in a dict.  
    '''
    return max(obj_dict.items(), key = lambda x: x[1])[0]




def sparse_mask(super_seg):
    ''' A mask to remove trivial values from intra-residue distances in a
        sparse matrix.
    '''
    # Load constant -- atomlabel...
    label_dict = constant_atomlabel()
    aa_dict    = constant_aminoacid_code()

    # Calculate the total length of distance matrix...
    len_list = [ len(label_dict[aa_dict[i]]) for i in super_seg ]
    len_dmat = np.sum( len_list )

    # Form a placeholder matrix with value one by default...
    dmask = np.zeros( (len_dmat, len_dmat))
    dmask[:] = 1

    # Assign zero to trivial values that only measure intra-residue distance...
    pos_current = 0
    for i in len_list:
        dmask[ pos_current : pos_current + i, pos_current : pos_current + i ] = 0.0
        pos_current += i

    return dmask
