#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pyrotein as pr
import givens as gv
from display import plot_dmat, plot_singular, plot_left_singular, plot_coeff
import multiprocessing as mp
from loaddata import load_xlsx, label_TMs
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
cterm = 322
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

if 0:
    # Visualize singular values...
    plot_singular(s, top = 10, log = True, index_from_zero = False)

if 0:
    labels_TM = label_TMs()
    for k, v in labels_TM.items(): labels_TM[k] = [ i * 4 for i in v ]

    # Visualize a u matrix...
    ## rank = 0
    ## plot_left_singular(u, rank, 
    ##                       length_mat = length_backbone, 
    ##                       guidelines = labels_TM,
    ##                       frac = 1.0,
    ##                       binning = 1)
    def plot_left_singualr_by_rank(rank):
        return plot_left_singular(u, rank, 
                                     length_mat = length_backbone, 
                                     guidelines = labels_TM,
                                     frac = 1.0,
                                     binning = 1,
                                     index_from_zero = False)
    with mp.Pool(num_job) as proc:
        proc.map( plot_left_singualr_by_rank, range(1, 11) )

if 1:
    # Create labels...
    entries = [ '-'.join(i[1:1+2]) for i in lines ]

    # Create color dictionary based on species...
    color_items = [ i[4] for i in lines ]

    # Define the color order according to a hypothesized reaction order...
    reaction_order = [ "11-cis", "11-cis detached", "9-cis", 
                        "batho", 
                        "lumi", 
                        "meta", 
                        "unobserved",
                        "All-trans detached",
                        "opsin" ]

    # Visualize...
    for j in range(1, 10):
        for i in range(j + 1,10):
            rank1, rank2 = j, i
            offset = "0.5,0"
            plot_coeff(c, rank1, rank2, entries = entries, 
                                        color_items = color_items, 
                                        color_order = reaction_order,
                                        label = True,
                                        ## xrange = (32.0, 32.5),
                                        ## yrange = (-0.5, 0.5),
                                        offset = offset, 
                                        rot = 0,
                                        height = 3,
                                        width = 3,
                                        index_from_zero = False)

if 0:
    # Create labels...
    point_labels = ['-'.join(i[1:1+2]) for i in lines]

    # Create color dictionary based on species...
    color_items = [ i[4] for i in lines ]

    # Define the color order according to a hypothesized reaction order...
    reaction_order = [ "11-cis", "11-cis detached", "9-cis", 
                        "batho", 
                        "lumi", 
                        "meta", 
                        "unobserved",
                        "All-trans detached",
                        "opsin" ]

    rank1, rank2 = 1, 4

    ## # Rotate points in rank1-rank2 plane by theta...
    ## theta = 90
    ## gv.givens_rotation(u, s, c, rank1, rank2, theta)

    offset = "0.5,0.0"
    plot_coeff(c, rank1, rank2, point_labels = point_labels, 
                                color_items = color_items, 
                                color_order = reaction_order,
                                ## xrange = (-0.8, -0.5),
                                ## yrange = (-0.2, 0.0),
                                offset = offset, 
                                rot = 0,
                                height = 3,
                                width = 3)

if 0:
    def plot_left_singualr_by_rank(rank):
        return plot_left_singular(u, rank, 
                                     length_mat = length_backbone,
                                     frac = 1.0,
                                     binning = 1), 
    with mp.Pool(num_job) as proc:
        proc.map( plot_left_singualr_by_rank, [rank1, rank2] )
