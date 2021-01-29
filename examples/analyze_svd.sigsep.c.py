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
    [2, 3,  28.6],
    [3, 4, -2.0],
    [4, 5,  -20],
    [5, 6,   40],
    [7, 8,   20],

    [2, 4,   0],
    [3, 5,  -4],
    [4, 6,  10],
    [5, 7,  20],
    [6, 8, -8.5],

    [2, 5,  20],

    [3, 8,  -4],
    [2, 7,  -20],
    [6, 7,   18],
    [8, 9,  26],
    [3, 7,  -20],
    [6, 9,  -30],
    [2, 6,   2.5],
]
for rank1, rank2, theta in rotations:
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)

# Create labels for u...
labels_TM = label_TMs()
for k, v in labels_TM.items(): labels_TM[k] = [ i * 4 for i in v ]

# Define the color order according to a hypothesized reaction order...
reaction_order = [ "11-cis", "11-cis detached", 
                    "unobserved", 
                    "9-cis", 
                    "13-cis",
                    "6mr",
                    "batho", 
                    "lumi", 
                    "All-trans detached",
                    "meta", 
                    "opsin" ]


# [[[ Customize ploting style ]]]
# Create style dictionary based on species...
coloritems = {}
coloritems["retinal"]  = {}
coloritems["methods"]  = {}
coloritems["comments"] = {}
coloritems["retinal"]  = select_items(lines, 4)
coloritems["methods"]  = select_items(lines, 9)
coloritems["comments"] = select_items(lines, 10)

# Define color by dividing the colorwheel...
reaction_colors_dict = cs.color_species(reaction_order, hexsym = "#")

# Assemble Gnuplot codes...
# retinal type
gps_dict = {}
for k, v in reaction_colors_dict.items():
    gps_dict[k] = { 
        "style" : "u 1:2 w p pt 7 ps 1 lc rgb '%s' title '%s'" % (v, k),
        "entry" : coloritems["retinal"][k]
    }

# Add support for EM structures
gps_dict["EM"] = {
    "style" : "u 1:2 w p pt 7 ps 0.4 lw 1 lc rgb 'black' title 'EM'",
    "entry" : coloritems["methods"]["Single particle EM"],
}

# Add support for stablized opsin
gps_dict["Bound to RS"] = {
    "style" : "u 1:2 w p pt 4 ps 0.4 lw 0.4 lc rgb 'black' title 'Bound to RS'",
    "entry" : coloritems["comments"]["Stabilizer"],
}

# Create labels...
entries = ['-'.join(i[1:1+2]) for i in lines]

cmds = [
       ## "unset xtics",
       ## "unset ytics",
       ## "unset xlabel",
       ## "unset ylabel",
       ## "unset border"
       f"set xzeroaxis",
       f"set yzeroaxis",
       ]

top = 10

if 1:
    if 0:
        # Visualize singular values...
        plot_singular(s, top = top, log = True, 
                      width = 5,
                      height = 5,
                      fontsize = 28,
                      linewidth = 2.0,
                      ticscale  = 2.0,
                      fl_path = f'eps_svd.rot',
                      index_from_zero = False)

    if 0:
        offset = "0.5,0.0"
        for j in range(1, top):
            for i in range(j + 1,top):
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
                           fl_path = f'eps_svd.rot',
                           fl_postfix = f'',
                           index_from_zero = False,
                           cmds = cmds)

    if 1:
        for rank in range(1, top):
            plot_blankcoeff(rank, width = 6, 
                                  height = 6, 
                                  fl_path = 'eps_svd.rot', 
                                  fl_postfix = f"", 
                                  index_from_zero = False)

    if 0:
        offset = "0.5,0.0"
        j = 2
        for i in range(j + 1,top):
            rank1, rank2 = j, i
            plot_coeff(c, rank2, rank1, entries = entries, 
                                        color_items = color_items, 
                                        color_order = reaction_order,
                                        color_saturation = 80,
                                        color_value = 100,
                                        label = True,
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
                                        fl_path = f'eps_svd.rot',
                                        fl_postfix = f'.v',
                                        index_from_zero = False,
                                        cmds = cmds)
