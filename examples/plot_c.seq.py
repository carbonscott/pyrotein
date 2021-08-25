#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pyrotein as pr
import givens as gv
from display import plot_dmat, plot_singular, plot_left_singular, \
                    plot_coeff, select_items, plot_blankcoeff
import multiprocessing as mp
from loaddata import load_xlsx
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

# Calculate the weighted coefficients...
# The following is a computational equivalence to c = RMS(u @ s, axis = 0) @ vh
c = np.matmul(np.diag(s), vh)
u_rms = 1.0 / np.sqrt(u.shape[0])
c     = c * u_rms

# Count #cpu for multiprocessing (optional)...
num_cpu = mp.cpu_count() // 2

# Define a series of rotation...
rotations = [
]
disp_index = -1    # 0-based Python convention
if len(rotations): rank1_last, rank2_last = rotations[disp_index][0:0 + 2]
for rank1, rank2, theta in rotations:
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)

# Define the color order according to a hypothesized reaction order...
species_list = [ 
                     "bovine", 
                     "squid", 
           ]

# [[[ Customize ploting style ]]]
# Create style dictionary based on species...
coloritems = {}
coloritems["state"]  = {}
coloritems["state"]  = select_items(lines, 3)

# Define color by dividing the colorwheel...
species_list1 = cs.color_species(species_list, s = 100, v = 100, hexsym = "#", b = 0, e = 90)
colors_dict = {**species_list1}

# Assemble Gnuplot codes...
# retinal type
gps_dict = {}
for k, v in colors_dict.items():
    gps_dict[k] = { 
        "style" : "u 1:2 w p pt 7 ps 2.0 lc rgb '%s' title '%s'" % (v, k),
        "entry" : coloritems["state"][k]
    }
cs.color_table(gps_dict)


# Create labels...
entries = ['-'.join(i[1:1+2]) for i in lines]

cmds = [
           f"set xzeroaxis",
           f"set yzeroaxis",
       ]

top = 10

if 1:
    if 1:
        # Visualize singular values...
        plot_singular(s, top = top, log = True, 
                      width = 6,
                      height = 5,
                      fontsize = 32,
                      linewidth = 1.0,
                      ticscale  = 2.0,
                      ## fl_path = f'svd.seq.c',
                      fl_path = f'figure.c',
                      index_from_zero = False)

    if 1:
        offset = "1.8,0.0"
        for j in range(1, top):
            ## for i in range(j + 1, top):
            for i in range(1, top):
                rank1, rank2 = j, i
                plot_coeff(c, rank1, rank2,
                           lbl = gps_dict,
                           ## label = True, 
                           label = False, 
                           labeltext = entries,
                           lbl_fontsize = 10,
                           xrange = (-1.5, 1.5),
                           yrange = (-1.5, 1.5),
                           offset = offset, 
                           rot = 0,
                           height = 6,
                           width = 6,
                           fontsize = 25,
                           pointsize = 2.0,
                           linewidth = 1.0,
                           fl_path = f'svd.seq.c',
                           ## fl_path = f'figure.c',
                           fl_postfix = f'',
                           index_from_zero = False,
                           cmds = cmds)


    if 0:
        offset = "1.8,0.0"
        for j in [rank1_last, rank2_last]:
            for i in range(1, top):
                rank1, rank2 = j, i
                plot_coeff(c, rank1, rank2,
                           plot_dict = gps_dict,
                           label = True, 
                           labeltext = entries,
                           lbl_fontsize = 10,
                           ## xrange = (0.25, 0.75),
                           ## yrange = (0.5, 0.9),
                           offset = offset, 
                           rot = 0,
                           height = 6,
                           width = 6,
                           fontsize = 25,
                           pointsize = 2.0,
                           linewidth = 1.0,
                           fl_path = f'svd.seq.c',
                           fl_postfix = f'',
                           index_from_zero = False,
                           cmds = cmds)


    if 0:
        for rank in range(1, top):
            plot_blankcoeff(rank, width = 6, 
                                  height = 6, 
                                  fl_path = 'svd.c', 
                                  fl_postfix = f"", 
                                  index_from_zero = False)


    if 0:
        offset = "1.8,0.0"
        cmds.append("set encoding utf8")
        cmds.append("set key top left")
        ## cmds.append("set key bot right")
        rank1, rank2 = 2, 5
        plot_coeff(c, rank1, rank2,
                   plot_dict = gps_dict,
                   ## label = True, 
                   label = False, 
                   labeltext = entries,
                   lbl_fontsize = 10,
                   ## xrange = (0.25, 0.75),
                   ## yrange = (0.5, 0.9),
                   offset = offset, 
                   rot = 0,
                   height = 6,
                   width = 6,
                   fontsize = 25,
                   pointsize = 2.0,
                   linewidth = 1.0,
                   ## fl_path = f'svd.c',
                   fl_path = f'figure.c',
                   fl_postfix = f'',
                   index_from_zero = False,
                   cmds = cmds)
