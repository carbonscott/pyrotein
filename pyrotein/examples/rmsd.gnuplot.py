#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pyrotein as pr
from display import plot_dmat

# Specify chains to process...
drc      = "pdb"
fl_chain = "chains.dat"
lines    = pr.utils.read_file(fl_chain)

# Define the backbone...
backbone = ["N", "CA", "C", "O"]

# Specify the range of atoms from rhodopsin...
nterm = 1
cterm = 348
len_backbone = (cterm - nterm + 1) * len(backbone)

# Initialize the matrix that stores accumulated coordinates...
dmats = np.zeros((len(lines), len_backbone, len_backbone))

# Accumulate coordinates...
for i_fl, (pdb, chain) in enumerate(lines):
    # Read coordinates from a PDB file...
    fl_pdb    = f"{pdb}.pdb"
    pdb_path  = os.path.join(drc, fl_pdb)
    atoms_pdb = pr.atom.read(pdb_path)

    # Create a lookup table for this pdb...
    atom_dict = pr.atom.create_lookup_table(atoms_pdb)

    # Obtain coordinates...
    xyzs = pr.atom.extract_backbone_xyz(atom_dict, chain, nterm, cterm)

    # Calculate individual distance matrix...
    dmat = pr.distance.calc_dmat(xyzs, xyzs)

    # Update the accumulated matrix...
    dmats[i_fl, :, :] = dmat[:, :]

# Calculate RMSD distance matrix...
rmsd_dmat = pr.distance.calc_rmsd_mats(dmats)

# Display the matrix...
rmsd_dmat_bin = pr.utils.bin_image(rmsd_dmat, bin = 8)
fl_rmsd_dmat  = "rmsd.eps"
pal = "set palette defined ( 0 'white', 0 'seagreen', 0.1 'white', 0.5 'blue', 1 'navy' )"
plot_dmat( rmsd_dmat_bin, 
           fl_rmsd_dmat, 
           lbl = [""] * len(rmsd_dmat_bin),
           lbl_font_size = 4,
           palette = pal,
           intst_max = 6,
           upper = 0.0,
           smooth = False,)
