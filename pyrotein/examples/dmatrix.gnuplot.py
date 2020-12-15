#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pyrotein as pr
from display import plot_dmat


# Read coordinates from a PDB file...
drc       =  "pdb"
pdb       =  "6cmo"
fl_pdb    = f"{pdb}.pdb"
chain     =  "R"
pdb_path  = os.path.join(drc, fl_pdb)
atoms_pdb = pr.atom.read(pdb_path)

# Create a lookup table for this pdb...
atom_dict = pr.atom.create_lookup_table(atoms_pdb)

# Define atoms used for distance matrix analysis...
peptide = ["N", "CA", "C", "O"]

# Specify the range of atoms from rhodopsin...
nterm = 1
cterm = 348

# Obtain the matrix of coordinates...
xyzs = pr.atom.extract_backbone_xyz(atom_dict, chain, nterm, cterm)

# Calculate the distance matrix...
dmat = pr.distance.calc_dmat(xyzs, xyzs)

# Export eps...
dmat_bin = np.array(pr.utils.bin_image(dmat, bin = 10))
fl_dmat = f"{pdb}.{chain}.dmat.eps"
pal = "set palette defined ( 0 'yellow', 0 'white', 0.5 'blue', 1 'navy' )"
plot_dmat(dmat_bin, fl_dmat, lbl = [""] * len(dmat_bin), palette = pal, smooth = True)
