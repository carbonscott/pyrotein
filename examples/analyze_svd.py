#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pyrotein as pr
import givens as gv
from display import plot_dmat, plot_singular, plot_left_singular, plot_coeff
import multiprocessing as mp
from loaddata import load_xlsx
import colorsimple as cs


def reverse_sign(u, vh, rank):
    u[:, rank]  = - u[:, rank]
    vh[rank, :] = -vh[rank, :]
    return None


# Specify chains to process...
fl_chain = "chains.comp.xlsx"
lines    = load_xlsx(fl_chain)

# Specify the range of atoms from rhodopsin...
nterm = 1
cterm = 348
backbone = ["N", "CA", "C", "O"]
length_backbone = (cterm - nterm + 1) * len(backbone)

# Load upstream data...
dmats = np.load("dmats.npy")
u     = np.load("u.npy")
s     = np.load("s.npy")
vh    = np.load("vh.npy")

# Allow positive value that fits the distance in nature (optinal)...
reverse_sign(u, vh, 0)

# Calculate the coefficients...
c = np.matmul(np.diag(s), vh)

# Standardize u and c and assign units...
u_ave = np.sqrt(u.shape[0])
c = c / u_ave

# Count #cpu for multiprocessing (optional)...
num_cpu = mp.cpu_count()
num_job = num_cpu - 1

if 1:
    # Visualize singular values...
    plot_singular(s, top = 10, log = True)

if 0:
    # Visualize a u matrix...
    ## for rank in range(1,10):
    ##     plot_left_singular(u, rank, 
    ##                           length_mat = length_backbone, 
    ##                           frac = 1.0,
    ##                           binning = 1)
    def plot_left_singualr_by_rank(rank):
        return plot_left_singular(u, rank, 
                                     length_mat = length_backbone,
                                     frac = 1.0,
                                     binning = 1)
    with mp.Pool(num_job) as proc:
        proc.map( plot_left_singualr_by_rank, range(0, 10) )

if 0:
    # Create labels...
    point_labels = [ '-'.join(i[1:1+2]) for i in lines ]

    # Create color dictionary based on species...
    color_items = [ i[4] for i in lines ]

    # Visualize...
    for j in range(0, 5):
        for i in range(j + 1,10):
            rank1, rank2 = j, i
            offset = "1.0,0"
            plot_coeff(c, rank1, rank2, point_labels = point_labels, 
                                        color_items = color_items, 
                                        ## xrange = (32.0, 32.5),
                                        ## yrange = (-0.5, 0.5),
                                        offset = offset, 
                                        rot = 0,
                                        height = 3,
                                        width = 3)

if 0:
    # Create labels...
    labels = ['_'.join(i[1:1+2]) for i in lines]

    rank1, rank2 = 0, 1

    # Rotate points in rank1-rank2 plane by theta...
    theta = 90
    gv.givens_rotation(u, s, c, rank1, rank2, theta)

    offset = "0.0,1.0"
    plot_coeff(c, rank1, rank2, labels = labels, offset = offset, rot = 90)

if 0:
    def plot_left_singualr_by_rank(rank):
        return plot_left_singular(u, rank, 
                                     length_mat = length_backbone,
                                     frac = 1.0,
                                     binning = 1), 
    with mp.Pool(num_job) as proc:
        proc.map( plot_left_singualr_by_rank, [rank1, rank2] )
