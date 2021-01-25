#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def rotation(v, theta):
    ''' Perform a spatial rotation.
        `theta` is in degree as an input.  
        `v` is a column vector with a dimension of 2 x 1.
        But a row vector works too.  
        It returns a vector with the same dimension.  
    '''
    theta = np.radians(theta)

    c = np.cos(theta)
    s = np.sin(theta)

    r = np.array( [ [c, -s], [s, c] ] )

    p = np.matmul(r, v)

    return p


def givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = True):
    ''' Perform a Givens rotation in the h-k plane around (0, 0)...
        `theta` is in degree.
    '''
    # Comply with the convention (1-based index)
    index1 = rank1 if index_from_zero else rank1 - 1
    index2 = rank2 if index_from_zero else rank2 - 1
    assert index1 >= 0, "Wrong value for rank1 (base-0 or base-1, is the input index correct?)"
    assert index2 >= 0, "Wrong value for rank2 (base-0 or base-1, is the input index correct?)"

    # Rotate c...
    c[(index1,index2), :] = rotation( c[(index1,index2), :], theta )

    # Find new s...
    u_ave = np.sqrt(u.shape[0])
    rescale_factor = u_ave * (u_ave - 1) / 2
    s[index1] = np.sqrt( np.sum( c[index1,:] * c[index1,:] ) * rescale_factor )
    s[index2] = np.sqrt( np.sum( c[index2,:] * c[index2,:] ) * rescale_factor )

    # Rotate u...
    u[:,(index1,index2)] = rotation( u[:,(index1,index2)].T, theta ).T

    return None
