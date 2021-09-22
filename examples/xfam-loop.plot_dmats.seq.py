#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
## sys.path.insert(0, "/home/scott/Dropbox/codes/pyrotein")
## sys.path.insert(0, "/Users/scott/Dropbox/codes/pyrotein")

import os
import numpy as np
import pyrotein as pr
from loaddata import load_gpcrdb_xlsx
from display import plot_dmat
import multiprocessing as mp

# Set job name...
job_name = "xfam-loop"

# [[[ OBTAIN THE CONSENSUS SEQUENCE ]]]
# Read the sequence alignment result...
fl_aln   = f"xfam.step4.psa.fil.fasta"
seq_dict = pr.fasta.read(fl_aln)

# Obtain the consensus sequence (super seq)...
tally_dict = pr.fasta.tally_resn_in_seqs(seq_dict)
super_seq  = pr.fasta.infer_super_seq(tally_dict)

# [[[ LOAD DATABASE ]]]
# Specify chains to process...
fl_chain = "gpcrdb.all.xlsx"
sheet    = f"xfam"
lines    = load_gpcrdb_xlsx(fl_chain, sheet = sheet, splitchain = True)
drc      = "pdb"

## lines = lines[:20]

# [[[ FIND SIZE OF DISTANCE MATRIX ]]]
# Get the sequence index (alignment) on the n-term side...
## nseqi = pr.fasta.get_lseqi(super_seq)
## cseqi = pr.fasta.get_rseqi(super_seq)

# The default ones are not good enough
nseqi = 637
cseqi = 1371

# Define atoms used for distance matrix analysis...
backbone = ["N", "CA", "C", "O"]
len_res = len(backbone)
len_seq = (cseqi - nseqi + 1) * len_res

# Accumulate distance matices as lower triangluar matrix...
len_lower_tri = (len_seq * len_seq - len_seq) // 2
dmats = np.zeros((len(lines), len_lower_tri))

# Build a mask to filter out loop segments...
mask_loop = np.zeros((len_seq, len("xyz")), dtype = np.int8)

# Use labels to create mask...
# Notice the +1 for right-inclusion is removed
labels = {'H8': [722, 734],
          'TM1': [0, 33],
          'TM2': [63, 95],
          'TM3': [158, 193],
          'TM4': [227, 253],
          'TM5': [316, 356],
          'TM6': [591, 630],
          'TM7': [680, 707]}
for k, v in labels.items(): labels[k] = [ (i) * len_res for i in v ]

# Create mask...
# +1 offset is to deal with numpy indexing which is not right-inclusive 
for _, (b, e) in labels.items(): mask_loop[b:e+1, :] = 1

# Preallocate memory for dmats...
dmats = np.zeros((len_seq, len_res, len_res))

drc_dmat = "xfam-loop.dmats"
pal = "set palette defined ( 0 '#F6FF9E', 0 'white', 0.5 'blue', 1 'navy' )"

if True:
    for i_fl, line in enumerate(lines[:]):
        # Unpack parameters
        uniprot = line[1].lower()
        spec    = line[5].lower()
        pdb     = line[7].lower()
        chain   = line[10].upper()

        if pdb != "7dd5": continue

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
        xyzs[mask_loop == 0] = np.nan

        # Calculate distance matrix...
        dmat = pr.distance.calc_dmat(xyzs, xyzs)

        fl_dmat = os.path.join(drc_dmat, f"{pdb}.{chain}.dmat")
        plot_dmat(dmat, fl_dmat, lbl = {}, palette = pal, NaN = 0)

        if pdb == "7dd5": break

if False:
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
