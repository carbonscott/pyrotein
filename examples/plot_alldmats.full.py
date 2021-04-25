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
import multiprocessing as mp

# [[[ OBTAIN THE CONSENSUS SEQUENCE ]]]
# Read the sequence alignment result...
# [WARNING] !!!sequence alignment is not trustworthy
fl_aln   = 'seq.align.fasta'
seq_dict = pr.fasta.read(fl_aln)

# Obtain the consensus sequence (super seq)...
tally_dict = pr.fasta.tally_resn_in_seqs(seq_dict)
super_seq  = pr.fasta.infer_super_seq(tally_dict)


# [[[ FIND SIZE OF DISTANCE MATRIX ]]]
# Get the sequence index (alignment) on the n-term side...
nseqi = pr.fasta.get_lseqi(super_seq)

# User defined range...
nterm, cterm = 1, 322
len_seg = cterm - nterm + 1
super_seg = super_seq[nseqi : nseqi + len_seg]

# Load constant -- atomlabel...
label_dict = pr.atom.constant_atomlabel()
aa_dict    = pr.atom.constant_aminoacid_code()

# [[[ ANALYZE PDB ENTRIES ]]]
# Specify chains to process...
fl_chain = "chains.comp.xlsx"
lines    = load_xlsx(fl_chain, sheet = "Sheet1")
drc      = "pdb"
drc_dmat = "dmats.full"

pal = "set palette defined ( 0 '#F6FF9E', 0 'white', 0.5 'blue', 1 'navy' )"
if True:
    for i_fl, line in enumerate(lines[-1:]):
        # Unpack parameters
        _, pdb, chain, _ = line[:4]

        # Read coordinates from a PDB file...
        fl_pdb    = f"{pdb}.pdb"
        pdb_path  = os.path.join(drc, fl_pdb)
        atoms_pdb = pr.atom.read(pdb_path)

        # Create a lookup table for this pdb...
        atom_dict = pr.atom.create_lookup_table(atoms_pdb)

        # Obtain the target protein by range...
        tar_seq = seq_dict[f"{pdb}_{chain}"]
        tar_seg = tar_seq[nseqi : nseqi + len_seg]

        # Standardize sidechain atoms...
        pr.atom.standardize_sidechain(atom_dict)

        # Obtain coordinates...
        xyzs = pr.atom.extract_xyz_by_seq(tar_seg, super_seg, atom_dict, chain, nterm, cterm)

        # Calculate distance matrix...
        dmat = pr.distance.calc_dmat(xyzs, xyzs)

        fl_dmat = os.path.join(drc_dmat, f"{pdb}.{chain}.dmat")
        ## plot_dmat(dmat, fl_dmat, lbl = {}, palette = pal, NaN = 0)
