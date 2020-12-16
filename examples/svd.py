#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pyrotein as pr
import os


# Specify chains to process...
fl_chain = "chains.comp.dat"
lines    = pr.utils.read_file(fl_chain)
drc      = "pdb"

# Define atoms used for distance matrix analysis...
backbone = ["N", "CA", "C", "O"]

# Specify the range of atoms from rhodopsin...
nterm = 1
cterm = 348
len_backbone = (cterm - nterm + 1) * len(backbone)

# Accumulate distance matices as lower triangluar matrix...
len_lower_tri = (len_backbone * len_backbone - len_backbone) // 2
dmats = np.zeros((len(lines), len_lower_tri))
for i_fl, (pdb, chain) in enumerate(lines):
    # Read coordinates from a PDB file...
    fl_pdb    = f"{pdb}.pdb"
    pdb_path  = os.path.join(drc, fl_pdb)
    atoms_pdb = pr.atom.read(pdb_path)

    # Create a lookup table for this pdb...
    atom_dict = pr.atom.create_lookup_table(atoms_pdb)

    # Obtain coordinates...
    xyzs = pr.atom.extract_backbone_xyz(atom_dict, chain, nterm, cterm)

    # Calculate distance matrix...
    dmat = pr.distance.calc_dmat(xyzs, xyzs)

    # Convert dmat into one-dimensional array and keep it in dmats...
    dmats[i_fl, :] = pr.utils.mat2tril(dmat, offset = -1)

# Replace np.nan with mean across samples...
pr.utils.fill_nan_with_mean(dmats.T, axis = 1)

# SVD...
# Column as example
# Row    as feature
u, s, vh = np.linalg.svd( dmats.T, full_matrices = False )
c = np.matmul(np.diag(s), vh)


# Export data for downstream analysis...
np.save("lines.npy" , lines)
np.save("dmats.npy" , dmats)
np.save("u.npy" , u)
np.save("s.npy" , s)
np.save("vh.npy", vh)
np.save("c.npy" , c)
