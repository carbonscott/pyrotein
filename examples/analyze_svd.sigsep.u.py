#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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


# Specify chains to process...
fl_chain = "chains.comp.xlsx"
lines    = load_xlsx(fl_chain)

# Specify the range of atoms from rhodopsin...
nterm = 1
cterm = 322    # It was 348
backbone = ["N", "CA", "C", "O"]
length_backbone = (cterm - nterm + 1) * len(backbone)

# Load upstream data...
dmats = np.load("dmats.npy")
u     = np.load("u.npy")
s     = np.load("s.npy")
vh    = np.load("vh.npy")

# Allow positive value that fits the distance in nature (optinal)...
reverse_sign(u, vh, 1, index_from_zero = False)
reverse_sign(u, vh, 2, index_from_zero = False)

# Calculate the coefficients...
c = np.matmul(np.diag(s), vh)

# Standardize u and c and assign units...
u_ave = np.sqrt(u.shape[0])
c = c / u_ave

# Count #cpu for multiprocessing (optional)...
num_cpu = mp.cpu_count() // 2

# Define a series of rotation...
rotations = [
    [2, 3, 17.5],
    [2, 5, 20.0],
    [3, 4, -0.0],
    [4, 5,  10],
    [5, 6, -32],
    [4, 6,  10],
    [8, 7,  -15],
    [8, 3,  6],
    [5, 7,  20],
    [2, 7, -20],
    [6, 7,  -15],
    [8, 9,  20],
]
disp_index = -6
if len(rotations): rank1_last, rank2_last = rotations[disp_index][0:0 + 2]
for rank1, rank2, theta in rotations:
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)

# Create labels...
entries = ['-'.join(i[1:1+2]) for i in lines]

# Create labels for u...
labels_TM = label_TMs()
for k, v in labels_TM.items(): labels_TM[k] = [ i * 4 for i in v ]

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

top = 10

if 1:
    # Visualize a u matrix...
    def plot_left_singualr_by_rank(rank):
        return plot_left_singular(u, rank, 
                                     length_mat = length_backbone, 
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
                                     fl_path = "eps_svd.rot",
                                     fl_postfix = f'',
                                     index_from_zero = False)
    if 0:
        num_job = np.min([top, num_cpu])
        with mp.Pool(num_job) as proc:
            proc.map( plot_left_singualr_by_rank, range(2, top) )
    if 1:
        num_job = 2
        with mp.Pool(num_job) as proc:
            proc.map( plot_left_singualr_by_rank, (rank1_last, rank2_last) )
