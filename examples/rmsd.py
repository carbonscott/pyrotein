#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## import sys
## sys.path.insert(0, "/home/scott/Dropbox/codes/pyrotein")
## sys.path.insert(0, "/Users/scott/Dropbox/codes/pyrotein")

import os
import numpy as np
import pyrotein as pr
from loaddata import load_xlsx, label_TMs
from display import plot_dmat


# Specify chains to process...
fl_chain = "chains.comp.xlsx"
lines    = load_xlsx(fl_chain)
drc      = "pdb"

# Define atoms used for distance matrix analysis...
peptide = ["N", "CA", "C", "O"]

# Specify the range of residues to analyze...
nterm = 1
cterm = 322
len_peptide = (cterm - nterm + 1) * len(peptide)

dmats = np.zeros((len(lines), len_peptide, len_peptide))
for i_fl, line in enumerate(lines):
    # Unpack parameters
    _, pdb, chain, species = line[:4]

    # Read coordinates from a PDB file...
    fl_pdb    = f"{pdb}.pdb"
    pdb_path  = os.path.join(drc, fl_pdb)
    atoms_pdb = pr.atom.read(pdb_path)

    # Create a lookup table for this pdb...
    atom_dict = pr.atom.create_lookup_table(atoms_pdb)

    # Obtain coordinates...
    xyzs = pr.atom.extract_xyz_by_atom(peptide, atom_dict, chain, nterm, cterm)

    # Calculate distance matrix...
    dmat = pr.distance.calc_dmat(xyzs, xyzs)
    dmats[i_fl, :, :] = dmat[:, :]

rmsd_dmat = pr.distance.calc_rmsd_mats(dmats)
np.save("rmsd_dmat.npy", rmsd_dmat)
