#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pyrotein as pr
from display import plot_dmat


# Read coordinates from a PDB file...
drc       =  "pdb_vr"
pdb       =  "4zwj"
fl_pdb    = f"{pdb}.pdb"
chain     =  "A"
pdb_path  = os.path.join(drc, fl_pdb)
atoms_pdb = pr.atom.read(pdb_path)

# Create a lookup table for this pdb...
atom_dict = pr.atom.create_lookup_table(atoms_pdb)

# Specify the range of atoms from rhodopsin...
nterm = 1
cterm = 348

# Obtain the coordinates...
xyzs = pr.atom.extract_backbone_xyz(atom_dict, chain, nterm, cterm)

# Calculate distance matrix...
dmat = pr.distance.calc_dmat(xyzs, xyzs)

# Export file...
## dmat_bin = dmat
## dmat_bin = np.array(pr.utils.bin_image(dmat, binning = 1))
dmat_bin = np.array(pr.utils.bin_image(dmat, binning = 2))
fl_dmat = f"{pdb}.{chain}.dmat"

# Colorscheme is inspired by from this publication (DOI: 10.1093/nar/gkw555) from Zhong Ren
pal = "set palette defined ( 0 '#F6FF9E', 0 'white', 0.5 'blue', 1 'navy' )"
lbl = [ f"" for i in range(len(dmat_bin)) ]
plot_dmat(dmat_bin, fl_dmat, lbl = lbl, palette = pal, smooth = True)
