#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pyrotein as pr
import givens as gv
from display import plot_dmat, plot_singular, plot_left_singular, \
                    plot_coeff, select_items, plot_blankcoeff
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
fl_path  = f"{job_name}.c"

# Specify chains to process...
fl_chain = "gpcrdb.all.xlsx"
sheet    = f"{job_name}"
lines    = load_gpcrdb_xlsx(fl_chain, sheet = sheet, splitchain = True)

# Load upstream data...
u     = np.load("u.seq.npy")
s     = np.load("s.seq.npy")
vh    = np.load("vh.seq.npy")
len_seq = np.load("len_seq.npy")

# Allow positive value that fits the distance in nature (optinal)...
reverse_sign(u, vh, 1, index_from_zero = False)
rev_list = [7]
for rev_i in rev_list: reverse_sign(u, vh, rev_i, index_from_zero = False)

# Calculate the weighted coefficients...
# The following is a computational equivalence to c = RMS(u @ s, axis = 0) @ vh
c = np.matmul(np.diag(s), vh)
u_rms = 1.0 / np.sqrt(u.shape[0])
c     = c * u_rms

# Define a series of rotation...
rotations = [
    [2,  3, 30.9-4.8],
    [2,  4, -8.1],
    [3,  4, 38.3],
    [3,  5, 14.8],
    [3,  6, -10.2],
    [4,  5, -24.4],
    [4,  6, 3.3],
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
if len(rotations): 
    rank1_last, rank2_last = rotations[disp_index][0:0 + 2]
for rank1, rank2, theta in rotations:
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)

## rank1_last, rank2_last = 2, 10


# Create style dictionary based on different coloring purposes...
# Register items from xlsx into coloritems
coloritems = {}

# Initialize an empty color table and plot statement...
colors_dict = {}
gps_dict    = {}

# Primarily hl by state...
hl_name = "state"
class_list = [ 
                 "A (Rhodopsin)", 
                 "B1 (Secretin)", 
                 "B2 (Adhesion)", 
                 "C (Glutamate)", 
                 "F (Frizzled)", 
              ]
class_list1 = cs.color_species(class_list, s = 100, v = 100, hexsym = "#", b = 0, e = 350)
class_list1["A (Rhodopsin)"] = "gray"
class_list1["B1 (Secretin)"] = "gray"
class_list1["B2 (Adhesion)"] = "gray"
class_list1["F (Frizzled)"] = "gray"

# Filter line by the activation state "Inactive"...
line_dict = {}
fil_name_dict = { "Inactive" : {"col" : 11, "shape" : 6, "lw" : 2}, 
                  "Active"   : {"col" : 11, "shape" : 7, "lw" : 1} }
is_exist_dict = {}
for fil_name, attr_dict in fil_name_dict.items():
    color_key = f"{hl_name} {fil_name}"
    col = attr_dict["col"]
    shape = attr_dict["shape"]
    lw = attr_dict["lw"]
    line_dict[color_key] = { i : line for i, line in enumerate(lines) if line[11] == fil_name }

    # States...
    coloritems[color_key] = select_items(line_dict[color_key], 4)
    colors_dict[color_key] = {**class_list1}
    for k, v in colors_dict[color_key].items():
        # Deal with wrongly-input receptor
        if not k in coloritems[color_key]:
            print(f"!!! {k} is not a valid recrod for coloring for style '{color_key}'.")
            continue

        # Set title...
        t = f"title '{fil_name} {k}'" if k == "C (Glutamate)" else f"title '{fil_name} Others'"
        if not t in is_exist_dict: is_exist_dict[t] = None
        else: t = "notitle"

        # Construct plot commands for each entry
        gps_dict[f"{color_key} {k}"] = { 
            "style" : f"u 1:2 w p pt {shape} ps 2.0 lw {lw} lc rgb '{v}' {t}",
            "entry" : coloritems[color_key][k]
        }

## # Activation...
## hl_name = "activation"
## hl_item = "Inactive"
## coloritems[hl_name]  = select_items(lines, 11)
## gps_dict[f"{hl_name}_{hl_item}"] = {
##     "style" : f"u 1:2 w p pt 6 ps 1.8 lw 2 lc rgb 'black' title '{hl_item}'",
##     "entry" : coloritems[hl_name][hl_item],
## }

cs.color_table(gps_dict)

# [[[ PLOT ]]]
# Create labels...
entry_dict = { f"{v[7]}_{v[10]}" : i for i, v in enumerate(lines) }

# Selectively plot entries...
pdb_list = ["5TE5_A", "1F88_A", "5TE3_A", "6K42_R", "6OY9_R", "6XBL_R",
            "6WW2_R", "5XEZ_B", "4L6R_A", "7CUM_B", "7C7Q_B", "4OR2_B",
            "7DD6_B", "6MXT_A", "3SN6_R", "4LDE_A", "6PS6_A", "2YCX_A",
            "2R4R_A", "6VCB_R", "6X1A_R", "6LML_R"]
entry_fil_dict = { k : entry_dict[k] for k in pdb_list if k in entry_dict }
## entry_fil_dict = entry_dict.copy()

cmds = [
           f"set xzeroaxis",
           f"set yzeroaxis",
       ]

top = 16

if 1:
    if 0:
        # Visualize singular values...
        plot_singular(s, top = top, log = True, 
                      width = 6,
                      height = 5,
                      fontsize = 32,
                      linewidth = 1.0,
                      ticscale  = 2.0,
                      fl_path = fl_path,
                      index_from_zero = False)

    if 0:
        offset = "1.8,0.0"
        start_index = 2
        for j in range(start_index, top):
            ## for i in range(start_index + 1, top):
            for i in range(1, top):
                rank1, rank2 = j, i
                plot_coeff(c, rank1, rank2,
                           lbl = gps_dict,
                           ## label = True, 
                           label = False, 
                           label_dict = entry_dict,
                           lbl_fontsize = 10,
                           xrange = (-2.0, 2.0),
                           yrange = (-2.0, 2.0),
                           offset = offset, 
                           rot = 0,
                           height = 6,
                           width = 6,
                           fontsize = 25,
                           pointsize = 2.0,
                           linewidth = 1.0,
                           fl_path = fl_path,
                           fl_postfix = f'',
                           index_from_zero = False,
                           cmds = cmds)

    if 0:
        cmds.append(f"set xtics 0.5")
        offset = "1.8,0.0"
        start_index = 3
        for j in range(start_index, top):
            rank1, rank2 = j, 1
            plot_coeff(c, rank1, rank2,
                       lbl = gps_dict,
                       label = False, 
                       label_dict = entry_dict,
                       lbl_fontsize = 10,
                       xrange = (-2.5, 2.5),
                       yrange = (23.6, 22.6),
                       offset = offset, 
                       rot = 0,
                       height = 6,
                       width = 6,
                       fontsize = 25,
                       pointsize = 2.0,
                       linewidth = 1.0,
                       fl_path = fl_path,
                       fl_postfix = f'',
                       index_from_zero = False,
                       cmds = cmds)

    if 1:
        offset = "1.8,0.0"
        cmds.append("set size ratio -1")
        rank1, rank2 = 7, 2
        plot_coeff(c, rank1, rank2,
                   lbl = gps_dict,
                   label = True, 
                   ## label_dict = entry_dict,
                   label_dict = entry_fil_dict,
                   lbl_fontsize = 10,
                   ## xrange = (-1.2, 1.8),
                   ## yrange = (-1.2, 1.8),
                   offset = offset, 
                   rot = 0,
                   height = 6,
                   width = 6,
                   fontsize = 25,
                   pointsize = 2.0,
                   linewidth = 1.0,
                   fl_path = fl_path,
                   fl_postfix = f'',
                   index_from_zero = False,
                   cmds = cmds)

    if 0:
        offset = "1.8,0.0"
        cmds.append("set encoding utf8")
        ## cmds.append("set key top left")
        start_index = 1
        rank1_adjust, rank2_adjust = rank1_last, rank2_last
        for i in [rank1_adjust, rank2_adjust]:
            for j in range(start_index, top):
                rank1, rank2 = i, j
                plot_coeff(c, rank1, rank2,
                           lbl = gps_dict,
                           ## label = True, 
                           label = False, 
                           label_dict = entry_fil_dict,
                           lbl_fontsize = 10,
                           xrange = (-3.0, 3.0),
                           yrange = (-3.0, 3.0),
                           offset = offset, 
                           rot = 0,
                           height = 6,
                           width = 6,
                           fontsize = 25,
                           pointsize = 2.0,
                           linewidth = 1.0,
                           fl_path = fl_path,
                           fl_postfix = f'',
                           index_from_zero = False,
                           cmds = cmds)

                rank1, rank2 = j, i
                plot_coeff(c, rank1, rank2,
                           lbl = gps_dict,
                           ## label = True, 
                           label = False, 
                           label_dict = entry_fil_dict,
                           lbl_fontsize = 10,
                           xrange = (-3.0, 3.0),
                           yrange = (-3.0, 3.0),
                           offset = offset, 
                           rot = 0,
                           height = 6,
                           width = 6,
                           fontsize = 25,
                           pointsize = 2.0,
                           linewidth = 1.0,
                           fl_path = fl_path,
                           fl_postfix = f'',
                           index_from_zero = False,
                           cmds = cmds)
