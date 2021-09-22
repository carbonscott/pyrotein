#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## import sys
## sys.path.insert(0, "/home/scott/Dropbox/codes/pyrotein")
## sys.path.insert(0, "/Users/scott/Dropbox/codes/pyrotein")

import os
import numpy as np
import pyrotein as pr
import givens as gv
from loaddata import load_gpcrdb_xlsx
from display import plot_dmat

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

# [ DEBUG ]
lines = lines[:4]

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
labels = {'H8': [722, 735],
          'TM1': [0, 34],
          'TM2': [63, 96],
          'TM3': [158, 194],
          'TM4': [227, 254],
          'TM5': [316, 357],
          'TM6': [591, 631],
          'TM7': [680, 708]}
for k, v in labels.items(): labels[k] = [ (i) * len_res for i in v ]

# Create mask...
for _, (b, e) in labels.items(): mask_loop[b:e, :] = 1.0

# Download the entrywise mean for dmats...
dmats_entrywise_mean = np.load(f"{job_name}.dmats.entrywise_mean.npy").reshape(-1)

# [[[ Obtain a structure ]]]
for i_fl, line in enumerate(lines):
    # Unpack parameters
    uniprot = line[1].lower()
    spec    = line[5].lower()
    pdb     = line[7].lower()
    chain   = line[10].upper()

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

    # Obtain xyzs
    entry = f"{pdb}_{chain}"
    print(entry)
    xyzs = pr.atom.extract_xyz_by_seq(backbone, chain_dict, tar_seq, nseqi, cseqi)
    xyzs[mask_loop == 0] = np.nan
    ## xyzs[mask_loop == 0] = 0.0

    # Calculate distance matrix...
    dmat = pr.distance.calc_dmat(xyzs, xyzs)
    dmat1d = pr.utils.mat2tril(dmat, offset = -1)

    # Patch nan based on column mean obtained prior to SVD...
    nan_pos = np.isnan(dmat1d)
    dmat1d[nan_pos] = dmats_entrywise_mean[nan_pos]
    dmats[i_fl, :] = dmat1d


# [[[ Load exisitng SVD heatmap U ]]]

def reverse_sign(u, vh, rank, index_from_zero = True):
    # Comply with the convention (1-based index)
    rank_in_data = rank if index_from_zero else rank - 1

    # Reverse sign...
    u[:, rank_in_data]  = - u[:, rank_in_data]
    vh[rank_in_data, :] = -vh[rank_in_data, :]

    return None

# Set job name...
aim      = "classA"
fl_path  = f"{job_name}.{aim}.c"


# Load upstream data...
b, e  = 0, 20 + 1
u       = np.load(f"{job_name}.u.seq.trunc.{b:02d}to{e:02d}.npy")
s       = np.load(f"{job_name}.s.seq.trunc.{b:02d}to{e:02d}.npy")
vh      = np.load(f"{job_name}.vh.seq.trunc.{b:02d}to{e:02d}.npy")
len_res = np.load(f"{job_name}.len_res.npy")
len_seq = np.load(f"{job_name}.len_seq.npy")

# Allow positive value that fits the distance in nature (optinal)...
reverse_sign(u, vh, 1, index_from_zero = False)
rev_list = [ ]
for rev_i in rev_list: reverse_sign(u, vh, rev_i, index_from_zero = False)

# Calculate the weighted coefficients...
# The following is a computational equivalence to c = RMS(u @ s, axis = 0) @ vh
c = np.matmul(np.diag(s), vh)
## u_rms = 1.0 / np.sqrt(u.shape[0])
u_rms = 1.0
u_rms = np.std(u[:, 1])
c     = c * u_rms

# Define a series of rotation...
rotations = [
   [ 2, 3, -45.3],
##   [ 5, 7,  19.9],
##   [ 2, 4,  13.7 + 3],
##   [ 2, 5, -11.0],
##   [ 4, 5, -14.6],
##   [ 4, 7, -39.5],
##   [ 5, 9,  11.2],
##   [ 6, 8,  59.9],
##   [ 3, 4,   5.1],
##   [ 3, 7, -13.5],
##   [ 3, 8,  21.6],
##   [ 5, 6,   4.2],
##   [ 5,10, - 9.1],
##   [ 4, 6,   3.4],
]
for i, (x, y, _) in enumerate(rotations):
    if x in rev_list: rotations[i][2] *= -1
    if y in rev_list: rotations[i][2] *= -1
disp_index = -1    # 0-based Python convention
if len(rotations):
    rank1_last, rank2_last = rotations[disp_index][0:0 + 2]
for rank1, rank2, theta in rotations:
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)



# [[[ Reconstruct coefficient ]]]

# The std of u...
## u_rms = 1.0 / np.sqrt(u.shape[0])

c_recstr = u.T @ dmats.T
c_recstr = c_recstr * u_rms

## # Treat nan as 0...
## dmat1d[np.isnan(dmat1d)] = 0.0
## 
## # Reconstruct...
## c_recstr = u.T @ dmat1d.reshape(-1,1) * u_rms
## vh_recstr = u.T @ dmat1d.reshape(-1,1)
