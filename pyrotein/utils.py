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




def read_file(file, str_to_num = False, num_type = float, labelcolumn = -1):
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
            if str_to_num: 
                words[labelcolumn + 1:] = [ num_type(word) for word in words[labelcolumn + 1:] ]
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




def sparse_mask(super_seg, offset = 1, val_offset = 0.0, val_init = 1.0):
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
    dmask[:] = val_init

    # Assign zero to trivial values that only measure intra-residue distance...
    len_resi = len(len_list)
    pos_x, pos_y = sum(len_list[ : offset]), 0
    for i, j in zip(len_list[ : len_resi - offset], len_list[ offset :]):
        dmask[ pos_x : pos_x + j, pos_y : pos_y + i ] = val_offset
        pos_x += j
        pos_y += i

    return dmask




def population_density(data, bin_cap = 100):
    ''' Return population density.
        bin_cap stands for bin capacity (number of items per bin).
    '''
    # Flatten data...
    data_flat = data.reshape(-1)

    # Sort data...
    data_sort = np.sort(data_flat)

    # Obtain the length of data...
    s, = data_sort.shape

    # Go through the array and figure out bin_val and bin_edge...
    bin_val  = []
    bin_edge = []
    bin_step = bin_cap
    for i in range(0, s, bin_cap):
        if i + bin_cap > s: bin_step = s - i
        data_seg = data_sort[i : i + bin_step]
        b, e     = data_seg[0], data_seg[-1]
        den      = bin_step / (e - b)
        bin_val.append(den)
        bin_edge.append(b)
    bin_edge.append( data_sort[-1] )

    return bin_val, bin_edge




def show_population_density(data, bin_cap, filename, 
                            rng       = [], 
                            width     = 3.5, 
                            height    = 2.62, 
                            fontsize  = 14, 
                            linewidth = 1.5,
                            xlabel    = 'Distance (\305)',
                            ylabel    = 'Population density (1/\305)',
                            linecolor = 'black',
                            cmds      = [],):
    data_val, data_rng = population_density(data, bin_cap = bin_cap)

    if len(rng) == 0: rng = data_rng[0], data_rng[-1]

    import GnuplotPy3
    gp = GnuplotPy3.GnuplotPy3()
    gp(f"set terminal postscript eps  size {width}, {height} \\")
    gp(f"                             enhanced color \\")
    gp(f"                             font 'Helvetica,{fontsize}' \\")
    gp(f"                             linewidth {linewidth}")
    gp(f"set output '{filename}.eps'")
    gp(f"set encoding utf8")
    gp(f"unset key")
    gp(f"set xlabel '{xlabel}'")
    gp(f"set ylabel '{ylabel}'")
    gp(f"set xrange [{rng[0]}:{rng[1]}]")

    for cmd in cmds:
        gp(cmd)

    gp("plot '-' using 1:2 with lines linewidth 1 linecolor rgb 'black'")

    for i in range(len(data_val)): 
        if data_rng[i] < rng[0]: continue
        if data_rng[i+1] > rng[1]: continue
        gp(f"{data_rng[i]} {data_val[i]}")  
        gp(f"{data_rng[i+1]} {data_val[i]}")  
    gp("e")

    gp("exit")

    return None




def label_dmat(super_seg, nterm, cterm):
    # Load constant -- atomlabel...
    label_dict = constant_atomlabel()
    aa_dict    = constant_aminoacid_code()

    # Go through residue and build a list of (resi, resn, atom)...
    label_list = []
    for seqi, resi in enumerate(range(nterm, cterm + 1)):
        aa        = super_seg[seqi]
        resn      = aa_dict[aa]
        atom_list = label_dict[resn]
        for atm in atom_list:
            label = f"{resi}.{resn}.{atm}"
            label_list.append(label)

    return label_list




def tally_list1d(int_list):
    int_dict = {}
    for i in int_list:
        if not i in int_dict: int_dict[i] = 1
        int_dict[i] += 1

    return int_dict




def chunker(seq, size = 60):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))




def sort_dict_by_key(_dict):
    ''' PDB doesn't sort resn by resi by default.  This function sorts resn by 
        resi in an ascending order.  
    '''
    return { k : v for k, v in sorted(_dict.items(), key = lambda x : x[0]) }




def sqeeze_seqi(lbl):
    ''' I know it's a terrible name.
        It turns things like this 
        [0, 33], [63, 95] -> 
        [0, 33], [34, 66]
    '''
    len_dict = { k : e - b for k, (b, e) in sorted(labels.items(), key = lambda x: [1][0]) }

    i_end = 0
    sqeezed_lbl = {}
    for k, l in len_dict.items():
        sqeezed_lbl[k] = [ i_end, i_end + l ]
        i_end += l+1

    return sqeezed_lbl
