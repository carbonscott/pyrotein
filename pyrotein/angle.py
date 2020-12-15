#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def calc_dihedral(b1, b2, b3):
    '''Refer to the formula on wikipedia: https://en.wikipedia.org/wiki/Dihedral_angle#cite_note-3
       for details.  
       It returns angles in a unit of degree.  

       b array has to be a 2D-array in order to support vectorization.

       >>> b = np.array([ [0.3, 0.5, 0.3] ])

       b array explained:
       - root axis     (axis 0): instances of b1 (or b2 or b3)
       - primary axisx (axis 1): XYZ coordinates
    '''
    b12     = np.cross(b1, b2, axis = 1)
    b23     = np.cross(b2, b3, axis = 1)
    b1213   = np.cross(b12, b23, axis = 1)
    b2_norm = b2 / np.linalg.norm(b2, axis = 1, keepdims = True)
    v1      = np.sum( b1213 * b2_norm, axis = 1 )
    v2      = np.sum( b12 * b23, axis =1 )

    return np.arctan2(v1, v2) / np.pi * 180.0
