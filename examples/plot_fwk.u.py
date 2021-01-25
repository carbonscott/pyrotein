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
    # Adjust the index base...
    index_base = 0 if index_from_zero else 1

    # Reverse sign...
    u[:, rank]  = - u[:, rank]
    vh[rank, :] = -vh[rank, :]

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
reverse_sign(u, vh, 4, index_from_zero = False)

# Calculate the coefficients...
c = np.matmul(np.diag(s), vh)

# Standardize u and c and assign units...
u_ave = np.sqrt(u.shape[0])
c = c / u_ave

# Count #cpu for multiprocessing (optional)...
num_cpu = mp.cpu_count()
num_job = num_cpu

# Draw guidelines (optional)...
lbls = {}
cmds_guideline_top = [""]
cmds_guideline_bottom = [""]
color_guideline = '#BBBBBB'

# Draw the rigid framework...
## fl_fwk = 'extracellular.fwk.dat'
fl_fwk = 'intracellular.fwk.dat'
fwk    = pr.utils.read_file(fl_fwk, numerical = True)
cmd_fwk = []
for i in range(len(fwk)):
    cmd = []
    b1, e1 = [ k * 4 for k in fwk[i] ]
    cmd.append(f"set object rectangle front from {b1},{e1} to {e1},{b1} fs empty border linecolor rgb 'black'")
    for j in range(i + 1, len(fwk)):
        b2, e2 = [ k * 4 for k in fwk[j] ]
        cmd.append(f"set object rectangle front from {b1},{b2} to {e1},{e2} fs empty linecolor rgb '{color_guideline}'")
    cmd_fwk.extend(cmd)
cmds_guideline_bottom.extend(cmd_fwk)
cmds_guideline_top.append(f"set key top left")

labels_TM = label_TMs()
for k, v in labels_TM.items(): labels_TM[k] = [ i * 4 for i in v ]

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
                                 cmds_top = cmds_guideline_top,
                                 cmds_bottom = cmds_guideline_bottom,
                                 index_from_zero = False)
plot_left_singualr_by_rank(2)
## num_job = 2
## with mp.Pool(num_job) as proc:
##     proc.map( plot_left_singualr_by_rank, [rank1, rank2] )


