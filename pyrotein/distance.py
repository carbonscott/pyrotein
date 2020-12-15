#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def calc_dmat(coor1, coor2):
    ''' Return a distance matrix.  
        coor1.shape: len(coor1), 3
        coor2.shape: len(coor2), 3
    '''
    # Convert them to np array...
    coor1_npary = np.array(coor1)
    coor2_npary = np.array(coor2)

    # Calculate the difference vector...
    delta_distance_vector = coor1_npary[:, None] - coor2_npary[None, :]

    # Calculate the distance matrix
    # np.sum is better suited than np.nansum for the purpose of distance matrix
    # For example, np.nansum([np.nan]) returns 0, whereas np.nan is a better outcome
    return np.sqrt( np.sum( delta_distance_vector * delta_distance_vector, axis = -1 ) )




def calc_rmsd_mats(mats):
    # Calculate the mean matrix...
    dmat_mean = np.nanmean( mats, axis = 0 )
    #                             ~~~~~~~~
    #                                 |
    # Collapse the 0th                |
    # dimension (Matrices)            |
    # with MEAN computation __________|

    # Difference matrix to obtain the deviation...
    # Make use of broadcasting w/o memory exploding
    dmats_diff = mats - dmat_mean[np.newaxis, :, :]
    #                             ~~~~~~~~~~  ~~~~
    #      New dimension ______________|        :
    # Original dimension .......................:

    # Obtain RMSD...
    dmat_rmsd = np.sqrt( np.nanmean( dmats_diff * dmats_diff, axis = 0 ) )

    return dmat_rmsd




def list_to_dmat(xyzs):
    """ Convert each list into a distance matrix.

        This function assumes more than one matrix is included.  For a single 
        matrix, use it as follows:

        list_to_dmat( [xyz] )
    """
    # Calculate distance matrix for each entry...
    dmats = np.zeros( (len(xyzs), len(xyzs[0]), len(xyzs[0])) )
    #                  ~~~~~~~~~  ~~~~~~~~~~~~  ~~~~~~~~~~~~
    # Entries _____________|            :             |
    # Rows .............................:             |
    # Cols ___________________________________________|

    for i in range(len(xyzs)):
        print(f"Computing matrix {i}...")
        dmat = calc_dmat( xyzs[i], xyzs[i] )
        dmats[i] = dmat.copy()
    return dmats




def tri_to_ary(xyzs):
    """ Convert the lower triangle of each matrix into a 1D array.

        This function assumes more than one matrix is included.  For a single 
        matrix, use it as follows:

        tri_to_ary( [xyz] )
    """
    # Size of lower triangle...
    len_mat  = len(xyzs[0])
    num_tril = ( len_mat * len_mat - len_mat ) // 2

    # Convert the lower triangle into 1D array...
    num_mat  = len(xyzs)
    data = np.zeros( (num_mat, num_tril) )
    for i in range(len(xyzs)):
        # Extract lower triangle w/o the main diag...
        tril = xyzs[i][np.tril_indices(len_mat,-1)]
        data[i] = tril
    return data
