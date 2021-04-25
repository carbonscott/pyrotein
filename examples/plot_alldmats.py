#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, "/home/scott/Dropbox/codes/pyrotein")
## sys.path.insert(0, "/Users/scott/Dropbox/codes/pyrotein")

import os
import numpy as np
import pyrotein as pr
from loaddata import load_xlsx, label_TMs
from display import plot_dmat
import multiprocessing as mp


# Specify chains to process...
fl_chain = "chains.comp.xlsx"
lines    = load_xlsx(fl_chain)
drc      = "pdb"

# Define atoms used for distance matrix analysis...
peptide = ["N", "CA", "C", "O"]

# Specify the range of atoms from adrenoceptor...
nterm = 1
cterm = 322

# The first element is to facilitate the indexing during assignment
len_segments = [ 0,
                 cterm - nterm + 1,
               ]
len_peptide = np.sum(len_segments) * len(peptide)

drc_dmat = "dmats"
pal = "set palette defined ( 0 '#F6FF9E', 0 'white', 0.5 'blue', 1 'navy' )"
if True:
    for i_fl, line in enumerate(lines[-1]):
        # Unpack parameters
        _, pdb, chain, species = line[:4]
        betatype = line[10]

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

        fl_dmat = os.path.join(drc_dmat, f"{pdb}.{chain}.dmat")
        plot_dmat(dmat, fl_dmat, lbl = {}, palette = pal, NaN = 0)


if False:
    def plotdmat(line):
        xyzs = np.zeros((len_peptide, 3))

        # Unpack parameters
        _, pdb, chain, species = line[:4]
        betatype = line[10]

        # Read coordinates from a PDB file...
        fl_pdb    = f"{pdb}.pdb"
        pdb_path  = os.path.join(drc, fl_pdb)
        atoms_pdb = pr.atom.read(pdb_path)

        # Create a lookup table for this pdb...
        atom_dict = pr.atom.create_lookup_table(atoms_pdb)

        # Obtain coordinates...
        xyzs = pr.atom.extract_xyz(peptide, atom_dict, chain, nterm, cterm)

        # Calculate distance matrix...
        dmat = pr.distance.calc_dmat(xyzs, xyzs)

        fl_dmat = os.path.join(drc_dmat, f"{pdb}.{chain}.dmat")
        plot_dmat(dmat, fl_dmat, lbl = {}, lbl_fontsize = 29, width = 10, height = 12, fontsize = 29, linewidth = 2.0, palette = pal, NaN = 0)

        return None


    num_job = 4
    with mp.Pool(num_job) as proc:
        proc.map( plotdmat, lines[:2] )
