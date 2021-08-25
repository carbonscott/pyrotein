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
fl_chain = "rhodopsin.db.xlsx"
lines    = load_xlsx(fl_chain, sheet = "total", splitchain = True)

# Load upstream data...
u     = np.load("u.seq.npy")
s     = np.load("s.seq.npy")
vh    = np.load("vh.seq.npy")
len_seq = np.load("len_seq.npy")

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
]
disp_index = -1    # 0-based Python convention
if len(rotations): rank1_last, rank2_last = rotations[disp_index][0:0 + 2]
for rank1, rank2, theta in rotations:
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)

# Create labels...
entries = ['-'.join(i[1:1+2]) for i in lines]

## # Create labels...
## # +1 to make it a right-end closed
## labels = {
##     "H1" : [  9,  32 + 1],
##     "H2" : [ 36,  62 + 1],
##     "H3" : [ 80, 101 + 1],
##     "H4" : [104, 128 + 1],
##     "H5" : [130, 161 + 1],
##     "H6" : [164, 192 + 1],
##     "H7" : [200, 225 + 1],
## }
## for k, v in labels.items(): labels[k] = [ (i - nterm) * len_res for i in v ]

cmds = [
       ## "unset xtics",
       ## "unset ytics",
       ## "unset xlabel",
       ## "unset ylabel",
       ## "unset border"
       ]

top = 1

if 1:
    # Visualize a u matrix...
    def plot_left_singualr_by_rank(rank):
        return plot_left_singular(u, rank, 
                                     length_mat = len_seq, 
                                     ## lbl = labels,
                                     width = 10,
                                     height = 12,
                                     fontsize = 29,
                                     lbl_fontsize = 29,
                                     linewidth = 2.0,
                                     frac = 1.0,
                                     binning = 1,
                                     intst_min = -0.01,
                                     intst_max =  0.01,
                                     fl_path = "svd.seq.u",
                                     fl_postfix = f'',
                                     index_from_zero = False)
    if 1:
        if __name__ == "__main__":
            for i in range(3):
                num_job = 4
                with mp.Pool(num_job) as proc:
                    proc.map( plot_left_singualr_by_rank, range(top + num_job * i, top + num_job * (i+1)) )
    if 0:
        num_job = 2
        with mp.Pool(num_job) as proc:
            proc.map( plot_left_singualr_by_rank, (rank1_last, rank2_last) )
