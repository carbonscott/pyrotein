#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def residue_complete(atom_dict, chain, backbone):
    ''' Select residues with a complete set of backbone atoms.
    '''

    resi_list = [ k for k, v in atom_dict[chain].items() \
                    if np.all([ i in v.keys() for i in backbone ]) ]

    return resi_list
