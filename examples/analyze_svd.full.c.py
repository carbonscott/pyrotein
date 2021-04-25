#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pyrotein as pr
import givens as gv
from display import plot_dmat, plot_singular, plot_left_singular, \
                    plot_coeff, select_items, plot_blankcoeff
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
lines    = load_xlsx(fl_chain, sheet = "Sheet1")

# Load upstream data...
dmats = np.load("dmats.full.npy")
u     = np.load("u.full.npy")
s     = np.load("s.full.npy")
vh    = np.load("vh.full.npy")

# Allow positive value that fits the distance in nature (optinal)...
reverse_sign(u, vh, 1, index_from_zero = False)
reverse_sign(u, vh, 2, index_from_zero = False)
reverse_sign(u, vh, 4, index_from_zero = False)
reverse_sign(u, vh, 6, index_from_zero = False)

# Calculate the coefficients...
c = np.matmul(np.diag(s), vh)

# Standardize u and c and assign units...
u_ave = np.sqrt(u.shape[0])
c = c / u_ave

# Count #cpu for multiprocessing (optional)...
num_cpu = mp.cpu_count() // 2

# Define a series of rotation...
rotations = [
    [3, 2, -25],
    [3, 4,  10],
    [5, 3,   7],
    [2, 4,   5],
    [5, 2, -10],
    [4, 6,  20],
    [8, 3,   8],
    [6, 5, -40],
    [5, 7,  18],
    [7, 6, -25],
    [7, 2,  20],
    [4, 3,   4],
    [3, 7,  28],
    [7, 6,  10],
    [9, 6,  30],
]
disp_index = -1    # 0-based Python convention
if len(rotations): rank1_last, rank2_last = rotations[disp_index][0:0 + 2]
for rank1, rank2, theta in rotations:
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)

# Create labels for u...
labels_TM = label_TMs()
for k, v in labels_TM.items(): labels_TM[k] = [ i * 4 for i in v ]

# Define the color order according to a hypothesized reaction order...
retinal_type = [ "11-cis", 
                 "11-cis detached", 
                 "unobserved", 
                 "9-cis", 
                 "6mr",
                 "batho", 
                 "lumi", 
                 "All-trans detached",
                 "13-cis",
                 "meta", 
                 "opsin", 
                 "drug binding"]


# [[[ Customize ploting style ]]]
# Create style dictionary based on species...
coloritems = {}
coloritems["retinal"]  = {}
coloritems["retinal"]  = select_items(lines, 4)

# Define color by dividing the colorwheel...
retinal_type_dict1 = cs.color_species(retinal_type[:8], s = 80, v = 100, hexsym = "#", 
                                                b = 20, e = 120)
retinal_type_dict2 = cs.color_species(retinal_type[-4:-1], s = 80, v = 100, hexsym = "#", 
                                                b = 180, e = 250)
retinal_type_dict3 = cs.color_species(retinal_type[-1:], s = 80, v = 100, hexsym = "#", 
                                                b = -60, e = 0)
retinal_colors_dict = {**retinal_type_dict1, **retinal_type_dict2, **retinal_type_dict3}

# Assemble Gnuplot codes...
# retinal type
gps_dict = {}
for k, v in retinal_colors_dict.items():
    gps_dict[k] = { 
        "style" : "u 1:2 w p pt 7 ps 2.0 lc rgb '%s' title '%s'" % (v, k),
        "entry" : coloritems["retinal"][k]
    }

# Add support for EM structures
coloritems["methods"]  = {}
coloritems["methods"]  = select_items(lines, 9)
gps_dict["EM"] = {
    "style" : "u 1:2 w p pt 6 ps 1.8 lw 2 lc rgb 'black' title 'EM'",
    "entry" : coloritems["methods"]["Single particle EM"],
}

# Arrestin binding
coloritems["Complex-C"] = {}
coloritems["Complex-C"] = select_items(lines, 7)
gps_dict["arrestin"] = {
    "style" : "u 1:2 w p pt 8 ps 1.5 lw 2 lc rgb 'black' title 'Arrestin binding'",
    "entry" : coloritems["Complex-C"]["S-arrestin:L374A/V375A/F376A"],
}

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
                      fl_path = f'svd.full.c',
                      index_from_zero = False)

    if 1:
        offset = "1.8,0.0"
        for j in range(1, top):
            ## for i in range(j + 1, top):
            for i in range(1, top):
                rank1, rank2 = j, i
                plot_coeff(c, rank1, rank2,
                           plot_dict = gps_dict,
                           label = True, 
                           ## label = False, 
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
                           fl_path = f'svd.full.c',
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
                           fl_path = f'svd.c',
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
