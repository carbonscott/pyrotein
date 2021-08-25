#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
## sys.path.insert(0, "/home/scott/Dropbox/codes/pyrotein")
## sys.path.insert(0, "/Users/scott/Dropbox/codes/pyrotein")

import os
import numpy as np
import pyrotein as pr
from loaddata import load_xlsx
from display import plot_dmat
import multiprocessing as mp

# [[[ OBTAIN THE CONSENSUS SEQUENCE ]]]
# Read the sequence alignment result...
# [WARNING] !!!sequence alignment is not trustworthy
fl_aln   = 'step3.msa.fasta'
seq_dict = pr.fasta.read(fl_aln)

# Obtain the consensus sequence (super seq)...
tally_dict = pr.fasta.tally_resn_in_seqs(seq_dict)
super_seq  = pr.fasta.infer_super_seq(tally_dict)

# [[[ FIND SIZE OF DISTANCE MATRIX ]]]
# Get the sequence index (alignment) on the n-term side...
nseqi = pr.fasta.get_lseqi(super_seq)
cseqi = pr.fasta.get_rseqi(super_seq)

# Specify chains to process...
fl_chain = "dopamine.db.xlsx"
lines    = load_xlsx(fl_chain, sheet = "total", splitchain = True)
drc      = "pdb"

# Define atoms used for distance matrix analysis...
peptide = ["N", "CA", "C", "O"]

# Define atoms used for distance matrix analysis...
backbone = ["N", "CA", "C", "O"]
len_res = len(backbone)
len_seq = (cseqi - nseqi + 1) * len_res

# Preallocate memory for dmats...
dmats = np.zeros((len_seq, len_res, len_res))

drc_dmat = "dmats.seq"
pal = "set palette defined ( 0 '#F6FF9E', 0 'white', 0.5 'blue', 1 'navy' )"

if False:
    for i_fl, line in enumerate(lines[:]):
        # Unpack parameters
        _, pdb, chain, nterm = line[:4]

        # Read coordinates from a PDB file...
        fl_pdb    = f"{pdb}.pdb"
        pdb_path  = os.path.join(drc, fl_pdb)
        atoms_pdb = pr.atom.read(pdb_path)

        # Create a lookup table for this pdb...
        atom_dict = pr.atom.create_lookup_table(atoms_pdb)

        # Obtain the chain to process...
        chain_dict = atom_dict[chain]

        # Obtain seq string for the current chain...
        tar_seq = seq_dict[f"{pdb.lower()}_{chain}"]

        # Obtain xyzs...
        entry = f"{pdb}_{chain}"
        xyzs = pr.atom.extract_xyz_by_seq(backbone, chain_dict, tar_seq, nseqi, cseqi)

        # Calculate distance matrix...
        dmat = pr.distance.calc_dmat(xyzs, xyzs)

        fl_dmat = os.path.join(drc_dmat, f"{pdb}.{chain}.dmat")
        plot_dmat(dmat, fl_dmat, lbl = {}, palette = pal, NaN = 0)


if True:
    def plotdmat(line):
        # Unpack parameters
        _, pdb, chain, lb_term, ub_term = line[:5]

        # Read coordinates from a PDB file...
        fl_pdb    = f"{pdb}.pdb"
        pdb_path  = os.path.join(drc, fl_pdb)
        atoms_pdb = pr.atom.read(pdb_path)

        # Create a lookup table for this pdb...
        atom_dict = pr.atom.create_lookup_table(atoms_pdb)

        # Obtain the chain to process...
        chain_dict = atom_dict[chain]

        # Obtain seq string for the current chain...
        tar_seq = seq_dict[f"{pdb.lower()}_{chain}"]

        # Obtain xyzs...
        entry = f"{pdb}_{chain}"
        xyzs = pr.atom.extract_xyz_by_seq(backbone, chain_dict, tar_seq, nseqi, cseqi)

        # Calculate distance matrix...
        dmat = pr.distance.calc_dmat(xyzs, xyzs)

        fl_dmat = os.path.join(drc_dmat, f"{pdb}.{chain}.dmat")
        plot_dmat(dmat, fl_dmat, lbl = {}, palette = pal, NaN = 0, linewidth = 1, intst_max = 100)

        return None


    if __name__ == "__main__":
        num_job = 4
        with mp.Pool(num_job) as proc:
            proc.map( plotdmat, lines[:] )
