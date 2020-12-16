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
