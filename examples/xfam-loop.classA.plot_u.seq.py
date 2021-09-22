#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pyrotein as pr
import givens as gv
from display import plot_dmat, plot_singular, plot_left_singular, plot_coeff
import multiprocessing as mp
from loaddata import load_gpcrdb_xlsx
import colorsimple as cs


def reverse_sign(u, vh, rank, index_from_zero = True):
    # Comply with the convention (1-based index)
    rank_in_data = rank if index_from_zero else rank - 1

    # Reverse sign...
    u[:, rank_in_data]  = - u[:, rank_in_data]
    vh[rank_in_data, :] = -vh[rank_in_data, :]

    return None


# Set job name...
job_name = "xfam-loop"
aim      = "classA"
fl_path  = f"{job_name}.{aim}.u"

# Specify chains to process...
fl_chain = "gpcrdb.all.xlsx"
sheet    = f"xfam"
lines    = load_gpcrdb_xlsx(fl_chain, sheet = sheet, splitchain = True)

# Load upstream data...
b, e    = 0, 20 + 1
u       = np.load(f"{job_name}.u.seq.trunc.{b:02d}to{e:02d}.npy")
s       = np.load(f"{job_name}.s.seq.trunc.{b:02d}to{e:02d}.npy")
vh      = np.load(f"{job_name}.vh.seq.trunc.{b:02d}to{e:02d}.npy")
len_res = np.load(f"{job_name}.len_res.npy")
len_seq = np.load(f"{job_name}.len_seq.npy")

# Allow positive value that fits the distance in nature (optinal)...
reverse_sign(u, vh, 1, index_from_zero = False)
rev_list = [ ]
for rev_i in rev_list: reverse_sign(u, vh, rev_i, index_from_zero = False)

# Calculate the coefficients...
c = np.matmul(np.diag(s), vh)

# Standardize u and c and assign units...
u_ave = np.sqrt(u.shape[0])
c = c / u_ave

# Define a series of rotation...
rotations = [
    [  2,  3, -48.8],
    [  2,  5, -14.5],
    [  3,  5, -43.4],
    [  5,  7,  43.4],
    [  8,  9,  43.4],
]
for i, (x, y, _) in enumerate(rotations):
    if x in rev_list: rotations[i][2] *= -1
    if y in rev_list: rotations[i][2] *= -1
disp_index = -1    # 0-based Python convention
if len(rotations): rank1_last, rank2_last = rotations[disp_index][0:0 + 2]
for rank1, rank2, theta in rotations:
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)

## rank1_last, rank2_last = 8, 11

# Create labels...
# +1 to make it a right-end closed
labels = {'H8': [722, 735],
          'TM1': [0, 34],
          'TM2': [63, 96],
          'TM3': [158, 194],
          'TM4': [227, 254],
          'TM5': [316, 357],
          'TM6': [591, 631],
          'TM7': [680, 708]}

# Increment the right index so that it is right-included
for k, v in labels.items(): labels[k][1] += 1

# Leave rooms for all atoms in a residue
for k, v in labels.items(): labels[k] = [ (i) * len_res for i in v ]


cmds = [
       ## "unset xtics",
       ## "unset ytics",
       ## "unset xlabel",
       ## "unset ylabel",
       ## "unset border"
       ]

top = 10


# Visualize a u matrix...
def plot_left_singualr_by_rank(rank):
    return plot_left_singular(u, rank, 
                                 length_mat = len_seq, 
                                 lbl = labels,

                                 ## width = 10,
                                 ## height = 12,
                                 ## fontsize = 29,
                                 ## lbl_fontsize = 29,

                                 width = 2.946,
                                 height = 3.535,
                                 fontsize = 8,
                                 lbl_fontsize = 8,

                                 lbl_linewidth = 1,
                                 curve_linewidth = 1,

                                 frac = 1.0,
                                 binning = 1,
                                 intst_min = -0.01,
                                 intst_max =  0.01,
                                 fl_path = fl_path,
                                 fl_postfix = f'',
                                 index_from_zero = False)

if 1:
    if __name__ == "__main__":
        num_job = 10
        if __name__ == "__main__":
            with mp.Pool(num_job) as proc:
                proc.map( plot_left_singualr_by_rank, range(1, 21) )
if 0:
    num_job = 2
    if __name__ == "__main__":
        with mp.Pool(num_job) as proc:
            proc.map( plot_left_singualr_by_rank, (rank1_last, rank2_last) )

if 0:
    for i in range(2, 16):
        plot_left_singualr_by_rank(i)

if 0:
    if __name__ == "__main__":
        num_job = 2
        with mp.Pool(num_job) as proc:
            proc.map( plot_left_singualr_by_rank, [rank1_last, rank2_last] )
