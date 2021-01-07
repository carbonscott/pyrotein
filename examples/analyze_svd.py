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

# Calculate the coefficients...
c = np.matmul(np.diag(s), vh)

# Standardize u and c and assign units...
u_ave = np.sqrt(u.shape[0])
c = c / u_ave

# Count #cpu for multiprocessing (optional)...
num_cpu = mp.cpu_count()
num_job = num_cpu

top = 10
if 0:
    # Visualize singular values...
    plot_singular(s, top = top, log = True, 
                  width = 5,
                  height = 5,
                  fontsize = 28,
                  linewidth = 2.0,
                  ticscale  = 2.0,
                  index_from_zero = False)

if 1:
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
                                     fl_path = "eps_svd",
                                     index_from_zero = False)
    ## num_job = 4
    ## with mp.Pool(num_job) as proc:
    ##     proc.map( plot_left_singualr_by_rank, range(1, top) )

    plot_left_singualr_by_rank(2)

if 1:
    # Create labels...
    entries = [ '-'.join(i[1:1+2]) for i in lines ]

    # Create color dictionary based on species...
    color_items = [ i[4] for i in lines ]

    # Define the color order according to a hypothesized reaction order...
    reaction_order = [ "11-cis", "11-cis detached", "9-cis", 
                        "batho", 
                        "lumi", 
                        "unobserved",
                        "All-trans detached",
                        "meta", 
                        "opsin", 
                     ]
    color_dict = cs.color_species(reaction_order)
    cs.color_table(color_dict)

    # Visualize...
    for j in range(1, top):
        for i in range(j + 1,top):
            rank1, rank2 = j, i
            offset = "0.5,0"
            plot_coeff(c, rank1, rank2, entries = entries, 
                                        color_items = color_items, 
                                        color_order = reaction_order,
                                        label = False,
                                        ## xrange = (32.0, 32.5),
                                        ## yrange = (-0.5, 0.5),
                                        lbl_fontsize = 10,
                                        offset = offset, 
                                        rot = 0,
                                        height = 6,
                                        width = 6,
                                        fontsize = 28,
                                        pointsize = 2.0,
                                        fl_path = 'eps_svd',
                                        index_from_zero = False)

            plot_coeff(c, rank1, rank2, entries = entries, 
                                        color_items = color_items, 
                                        color_order = reaction_order,
                                        label = True,
                                        ## xrange = (32.0, 32.5),
                                        ## yrange = (-0.5, 0.5),
                                        lbl_fontsize = 10,
                                        offset = offset, 
                                        rot = 0,
                                        height = 6,
                                        width = 6,
                                        fontsize = 28,
                                        pointsize = 2.0,
                                        fl_path = 'eps_svd',
                                        index_from_zero = False)
# Zoom or Rotate
if 0:
    # Create labels...
    entries = ['-'.join(i[1:1+2]) for i in lines]

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

    rank1, rank2 = 3, 4

    # Rotate points in rank1-rank2 plane by theta...
    theta = 335
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)

    cmds = [
           ## "unset xtics",
           ## "unset ytics",
           ## "unset xlabel",
           ## "unset ylabel",
           ## "unset border"
           ]

    offset = "0.5,0.0"
    plot_coeff(c, rank1, rank2, entries = entries, 
                                color_items = color_items, 
                                color_order = reaction_order,
                                label = True,
                                ## xrange = (0.25, 0.75),
                                ## yrange = (0.5, 0.9),
                                offset = offset, 
                                rot = 0,
                                height = 3,
                                width = 3,
                                linewidth = 1.0,
                                index_from_zero = False,
                                cmds = cmds)

# Rotate u
if 0:
    labels_TM = label_TMs()
    for k, v in labels_TM.items(): labels_TM[k] = [ i * 4 for i in v ]

    # Visualize a u matrix...
    def plot_left_singualr_by_rank(rank):
        return plot_left_singular(u, rank, 
                                     length_mat = length_backbone, 
                                     guidelines = labels_TM,
                                     fontsize = 14,
                                     lbl_fontsize = 28,
                                     linewidth = 2.0,
                                     frac = 1.0,
                                     binning = 1,
                                     index_from_zero = False)
    num_job = 2
    with mp.Pool(num_job) as proc:
        proc.map( plot_left_singualr_by_rank, [rank1, rank2] )


# Rotate c
if 0:
    # Create labels...
    entries = ['-'.join(i[1:1+2]) for i in lines]

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

    rank1, rank2 = 3, 4

    cmds = [
           ## "unset xtics",
           ## "unset ytics",
           ## "unset xlabel",
           ## "unset ylabel",
           ## "unset border"
           ]
    offset = "0.5,0.0"

    # Rotate points in rank1-rank2 plane by theta...
    ## theta = 10
    delta_theta = 5
    for theta in range(0,360,1):
        gv.givens_rotation(u, s, c, rank1, rank2, delta_theta, index_from_zero = False)

        plot_coeff(c, rank1, rank2, entries = entries, 
                                    color_items = color_items, 
                                    color_order = reaction_order,
                                    label = True,
                                    ## xrange = (0.25, 0.75),
                                    ## yrange = (0.5, 0.9),
                                    offset = offset, 
                                    rot = 0,
                                    height = 3,
                                    width = 3,
                                    linewidth = 1.0,
                                    index_from_zero = False,
                                    cmds = cmds)

        print(f"Current rotation: {delta_theta * theta:03.2f} -- ", end = "")
        opt = input("[e]xit or any key to continue:").lower()
        if opt == "e": break

if 0:
    labels_TM = label_TMs()
    for k, v in labels_TM.items(): labels_TM[k] = [ i * 4 for i in v ]

    # Visualize a u matrix...
    def plot_left_singualr_by_rank(rank):
        return plot_left_singular(u, rank, 
                                     length_mat = length_backbone, 
                                     guidelines = labels_TM,
                                     fontsize = 14,
                                     lbl_fontsize = 14,
                                     frac = 1.0,
                                     binning = 1,
                                     index_from_zero = False)
    with mp.Pool(num_job) as proc:
        proc.map( plot_left_singualr_by_rank, [rank1, rank2] )

# Partial
if 0:
    labels_TM = label_TMs()
    for k, v in labels_TM.items(): labels_TM[k] = [ i * 4 for i in v ]

    vrange = [199, 277]
    labels_vrange = {}
    for i in range(vrange[0], vrange[1]):
        labels_vrange[f"{i}"] = [ i * 4, i * 4 + 4]
    vrange = [ len(backbone) * i for i in vrange ]

    rank = 2
    plot_left_singular(u, rank, 
                          length_mat = length_backbone, 
                          ## guidelines = labels_vrange,
                          guidelines = labels_TM,
                          showguidelines = True,
                          width = 6,
                          height= 7,
                          fontsize = 14,
                          lbl_fontsize = 14,
                          frac = 1.0,
                          binning = 1,
                          ## vrange  = vrange,
                          index_from_zero = False)

