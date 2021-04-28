#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## import sys
## sys.path.insert(0, "/home/scott/Dropbox/codes/pyrotein")

import numpy as np
import pyrotein as pr
import givens as gv
from display import plot_dmat, plot_singular, plot_left_singular, plot_coeff
import multiprocessing as mp
from loaddata import load_xlsx, label_TMs
import colorsimple as cs


def reverse_sign(u, vh, rank, index_from_zero = True):
    # Comply with the convention (1-based index)
    rank_in_data = rank if index_from_zero else rank - 1

    # Reverse sign...
    u[:, rank_in_data]  = - u[:, rank_in_data]
    vh[rank_in_data, :] = -vh[rank_in_data, :]

    return None

length_mat = np.load("len_dmat.npy")

# Specify chains to process...
fl_chain = "chains.comp.xlsx"
lines    = load_xlsx(fl_chain, sheet = "Sheet1")

# Load upstream data...
dmats_sub = np.load("dmats.contact.npy")
u     = np.load("u.contact.npy")
s     = np.load("s.contact.npy")
vh    = np.load("vh.contact.npy")
in_range_index = np.load("index.contact.npy")
len_lower_tri = np.load("len_lower_tri.npy")

# Allow positive value that fits the distance in nature (optinal)...
reverse_sign(u, vh, 1, index_from_zero = False)
## reverse_sign(u, vh, 2, index_from_zero = False)
## reverse_sign(u, vh, 4, index_from_zero = False)
## reverse_sign(u, vh, 6, index_from_zero = False)

# Calculate the coefficients...
c = np.matmul(np.diag(s), vh)

# Standardize u and c and assign units...
u_ave = np.sqrt(u.shape[0])
c = c / u_ave

# Count #cpu for multiprocessing (optional)...
num_cpu = mp.cpu_count() // 2

# Define a series of rotation...
rotations = [
    [ 2, 3, 20 ],
    [ 3, 4, -3 ],
    [ 4, 6, 35 ],
    [ 4, 5, 10 ],
    [ 5, 7, 40 ],
]
disp_index = -1    # 0-based Python convention
if len(rotations): rank1_last, rank2_last = rotations[disp_index][0:0 + 2]
for rank1, rank2, theta in rotations:
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)

# Create labels...
entries = ['-'.join(i[1:1+2]) for i in lines]

# [[[ OBTAIN THE CONSENSUS SEQUENCE ]]]
# Read the sequence alignment result...
# [WARNING] !!!sequence alignment is not trustworthy, need to check manually
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

# Create labels for u...
labels_TM = label_TMs()
for k, v in labels_TM.items(): 
    labels_TM[k] = [ pr.fasta.resi_to_seqi(i, super_seg, nterm) for i in v ]

# Create color dictionary based on species...
color_items = [ i[4] for i in lines ]

# Define the color order according to a hypothesized reaction order...
reaction_order = [ "11-cis", "11-cis detached", "9-cis", 
                    "batho", 
                    "lumi", 
                    "unobserved",
                    "All-trans detached",
                    "meta", 
                    "opsin" ]

cmds = [
       ## "unset xtics",
       ## "unset ytics",
       ## "unset xlabel",
       ## "unset ylabel",
       ## "unset border"
       ]

# Restore sparse matrix
u_full = np.zeros((len(lines), len_lower_tri))
u_full[:] = np.nan
u_full[:, in_range_index] = u.T

top = 10

if 1:
    # Visualize a u matrix...
    def plot_left_singualr_by_rank(rank):
        ## return plot_left_singular(u, rank, 
        return plot_left_singular(u_full.T, rank, 
                                     length_mat = length_mat, 
                                     guidelines = labels_TM,
                                     width = 10,
                                     height = 12,
                                     fontsize = 29,
                                     lbl_fontsize = 29,
                                     linewidth = 2.0,
                                     frac = 1.0,
                                     binning = 1,
                                     intst_min = -0.01,
                                     intst_max =  0.01,
                                     fl_path = "svd.contact.u",
                                     fl_postfix = f'',
                                     quick = True,
                                     index_from_zero = False)
    if 1:
        num_job = np.min([top, num_cpu])
        with mp.Pool(num_job) as proc:
            proc.map( plot_left_singualr_by_rank, range(1, top + 1) )
    if 0:
        num_job = 2
        with mp.Pool(num_job) as proc:
            proc.map( plot_left_singualr_by_rank, (rank1_last, rank2_last) )
