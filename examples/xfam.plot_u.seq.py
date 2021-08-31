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
job_name = "xfam"
fl_path  = f"{job_name}.u"

# Specify chains to process...
fl_chain = "gpcrdb.all.xlsx"
sheet    = f"{job_name}"
lines    = load_gpcrdb_xlsx(fl_chain, sheet = sheet, splitchain = True)

# Load upstream data...
u     = np.load("u.seq.npy")
s     = np.load("s.seq.npy")
vh    = np.load("vh.seq.npy")
len_res = np.load("len_res.npy")
len_seq = np.load("len_seq.npy")

# Allow positive value that fits the distance in nature (optinal)...
reverse_sign(u, vh, 1, index_from_zero = False)
rev_list = [7]
for rev_i in rev_list: reverse_sign(u, vh, rev_i, index_from_zero = False)

# Calculate the coefficients...
c = np.matmul(np.diag(s), vh)

# Standardize u and c and assign units...
u_ave = np.sqrt(u.shape[0])
c = c / u_ave

# Define a series of rotation...
rotations = [
    [2, 3, 30.9-4.8],
    [2, 4, -8.1],
    [3, 4, 38.3],
    [3, 5, 14.8],
    [3, 6, -10.2],
    [4, 5, -24.4],
    [4, 6, 3.3],
    [4, 11, -41.8],
    [6,  8, -41.0],
    [6,  7,  25.1],
    [8,  2, -10.2],
    [8,  7, 18.7],
    [8,  9,-11.1],
    [9,  6,-33.5],
    [9,  7,-14.3],
    [9, 10,-22.9],
    [7, 11,-18.1],
    [7,  5,-22.8],
    [7,  4,-11.4],
    [2,  7,  5.6],
]
for i, (x, y, _) in enumerate(rotations):
    if x in rev_list: rotations[i][2] *= -1
    if y in rev_list: rotations[i][2] *= -1
disp_index = -1    # 0-based Python convention
if len(rotations): rank1_last, rank2_last = rotations[disp_index][0:0 + 2]
for rank1, rank2, theta in rotations:
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)

## rank1_last, rank2_last = 6, 10

# Create labels...
# +1 to make it a right-end closed
labels = {'H8': [620, 633],
          'TM1': [0, 33],
          'TM2': [47, 79],
          'TM3': [115, 151],
          'TM4': [182, 209],
          'TM5': [277, 317],
          'TM6': [503, 542],
          'TM7': [584, 612]}
for k, v in labels.items(): labels[k] = [ (i) * len_res for i in v ]

cmds = [
       ## "unset xtics",
       ## "unset ytics",
       ## "unset xlabel",
       ## "unset ylabel",
       ## "unset border"
       ]

top = 1


# Visualize a u matrix...
def plot_left_singualr_by_rank(rank):
    return plot_left_singular(u, rank, 
                                 length_mat = len_seq, 
                                 lbl = labels,
                                 width = 10,
                                 height = 12,
                                 fontsize = 29,
                                 lbl_fontsize = 29,
                                 frac = 1.0,
                                 binning = 1,
                                 intst_min = -0.01,
                                 intst_max =  0.01,
                                 fl_path = fl_path,
                                 fl_postfix = f'',
                                 index_from_zero = False)

if 0:
    if __name__ == "__main__":
        for i in range(16):
            num_job = 1
            with mp.Pool(num_job) as proc:
                proc.map( plot_left_singualr_by_rank, range(top + num_job * i, top + num_job * (i+1)) )
if 0:
    num_job = 2
    with mp.Pool(num_job) as proc:
        proc.map( plot_left_singualr_by_rank, (rank1_last, rank2_last) )

if 0:
    for i in range(2, 16):
        plot_left_singualr_by_rank(i)

if 1:
    rank1, rank2 = rank1_last, rank2_last
    plot_left_singualr_by_rank(rank1)
    plot_left_singualr_by_rank(rank2)
