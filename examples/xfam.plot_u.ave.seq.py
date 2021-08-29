#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pyrotein as pr
import givens as gv
from display import plot_u_ave
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
fl_path  = f"{job_name}.ave"

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
## reverse_sign(u, vh, 2, index_from_zero = False)
## reverse_sign(u, vh, 4, index_from_zero = False)
## reverse_sign(u, vh, 6, index_from_zero = False)

# Calculate the coefficients...
c = np.matmul(np.diag(s), vh)

# Standardize u and c and assign units...
u_ave = np.sqrt(u.shape[0])
c = c / u_ave

# Define a series of rotation...
rotations = [
]
disp_index = -1    # 0-based Python convention
if len(rotations): rank1_last, rank2_last = rotations[disp_index][0:0 + 2]
for rank1, rank2, theta in rotations:
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)

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

if 1:
    # Visualize a u matrix...
    def plot_u_ave_by_rank(rank):
        return plot_u_ave( u,
                            rank,
                            ## lbl = labels,
                            xrange = (750, 1000),
                            fl_postfix = "zoom",
                            yrange = (-0.0016, 0.0016),
                            length_mat = len_seq,
                            fl_path = fl_path,
                            index_from_zero = False)

    if 1:
        plot_u_ave_by_rank(2)

    if 0:
        for i in range(2, 16):
            plot_u_ave_by_rank(i)
    if 0:
        num_job = 2
        with mp.Pool(num_job) as proc:
            proc.map( plot_left_singualr_by_rank, (rank1_last, rank2_last) )
