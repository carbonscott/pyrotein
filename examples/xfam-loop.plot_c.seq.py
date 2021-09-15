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
from itertools import product, permutations


def reverse_sign(u, vh, rank, index_from_zero = True):
    # Comply with the convention (1-based index)
    rank_in_data = rank if index_from_zero else rank - 1

    # Reverse sign...
    u[:, rank_in_data]  = - u[:, rank_in_data]
    vh[rank_in_data, :] = -vh[rank_in_data, :]

    return None


# Set job name...
job_name = "xfam-loop"
fl_path  = f"{job_name}.c"

# Specify chains to process...
fl_chain = "gpcrdb.all.xlsx"
sheet    = f"xfam"
lines    = load_gpcrdb_xlsx(fl_chain, sheet = sheet, splitchain = True)

# Load upstream data...
b, e  = 0, 40 + 1
u       = np.load(f"{job_name}.u.seq.trunc.{b:02d}to{e:02d}.npy")
s       = np.load(f"{job_name}.s.seq.trunc.{b:02d}to{e:02d}.npy")
vh      = np.load(f"{job_name}.vh.seq.trunc.{b:02d}to{e:02d}.npy")
len_res = np.load(f"{job_name}.len_res.npy")
len_seq = np.load(f"{job_name}.len_seq.npy")

# Allow positive value that fits the distance in nature (optinal)...
reverse_sign(u, vh, 1, index_from_zero = False)
rev_list = [4]
for rev_i in rev_list: reverse_sign(u, vh, rev_i, index_from_zero = False)

# Calculate the weighted coefficients...
# The following is a computational equivalence to c = RMS(u @ s, axis = 0) @ vh
c = np.matmul(np.diag(s), vh)
u_rms = 1.0 / np.sqrt(u.shape[0])
c     = c * u_rms

# Define a series of rotation...
rotations = [
    [ 3, 4, -32.5],
    [ 2, 4,  51.8],
    [ 2, 5, -14.4],
]
for i, (x, y, _) in enumerate(rotations):
    if x in rev_list: rotations[i][2] *= -1
    if y in rev_list: rotations[i][2] *= -1
disp_index = -1    # 0-based Python convention
if len(rotations):
    rank1_last, rank2_last = rotations[disp_index][0:0 + 2]
for rank1, rank2, theta in rotations:
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)

## rank1_last, rank2_last = 3, 5


# Create style dictionary based on different coloring purposes...
# Register items from xlsx into coloritems
coloritems = {}

# Initialize an empty color table and plot statement...
colors_dict = {}
gps_dict    = {}

# Primarily hl by Classes...
class_list = [
                 ## "A (Rhodopsin)",
                 "B1 (Secretin)",
                 "B2 (Adhesion)",
                 "C (Glutamate)",
                 "F (Frizzled)",
                 ## "D1 (Ste2-like fungal pheromone)",
              ]
class_list1 = {}
## class_list1["A (Rhodopsin)"] = "#FFB8B8"    # Low key red
## class_list1["B1 (Secretin)"] = "#00A600"
## class_list1["C (Glutamate)"] = "#0080FF"
## class_list1["B2 (Adhesion)"] = "gray"
## class_list1["F (Frizzled)"] = "#AB00FF"

class_list1["B1 (Secretin)"] = "#E6E6E6"
class_list1["C (Glutamate)"] = "#E6E6E6"
class_list1["B2 (Adhesion)"] = "#E6E6E6"
class_list1["F (Frizzled)"]  = "#E6E6E6"

# Filter line by the activation state "Inactive"...
line_dict = {}
state_dict = {
                  "Inactive" : {"col" : 11, "shape" : 6, "lw" : 1},
                  "Intermediate" : {"col" : 11, "shape" : 7, "lw" : 1, },
                  "Active"   : {"col" : 11, "shape" : 2, "lw" : 1}
                }
is_exist_dict = {}
hl_name = "Classes"
for fil_name, attr_dict in state_dict.items():
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
            print(f"!!! {k} is not a valid record for coloring for style '{color_key}'.")
            continue

        # Set title...
        # Title will affect the making of key
        ## t = f"title '{fil_name} {k}'" if k == "C (Glutamate)" or k == "F (Frizzled)" else f"title '{fil_name} Others'"

        ## t = f"title '{fil_name} {k}'"

        t = f"title '{k}'"

        if not t in is_exist_dict: is_exist_dict[t] = None
        else: t = "notitle"

        # Construct plot commands for each entry
        gps_dict[f"{color_key} {k}"] = {
            "style" : f"u 1:2 w p pt {shape} ps 1.0 lw {lw} lc rgb '{v}' {t}",
            "entry" : coloritems[color_key][k]
        }

## gps_dict = {}

classAfamily_tier_dict = {
    "tier 1": [
        "Adrenoceptors",
        "Adenosine",
        "Opsins",
    ],

    "tier 2" : [
        "5-Hydroxytryptamine",
        "Orexin",
        "Acetylcholine (muscarinic)",
        "Dopamine",
        "Chemokine",
        "Neurotensin",
        "Opioid",
    ],

    "tier 3" : [
        "Prostanoid",
        "Angiotensin",
        "Cannabinoid",
        "Melatonin",
    ],

    "tier 4" : [
        "Ghrelin",
        "Platelet-activating factor",
        "Lysophospholipid (S1P)",
        "Apelin",
        "Complement peptide",
        "Cholecystokinin",
        "Formylpeptide",
        "Histamine",
        "Melanocortin",
        "Bile acid",
        "Succinate",
        "Vasopressin and oxytocin",
        "Neuropeptide Y",
        "Free fatty acid",
        "A orphans",
        "Lysophospholipid (LPA)",
        "Proteinase-activated",
        "P2Y",
        "Tachykinin",
        "Leukotriene",
        "Endothelin",
    ],
}
classAfamily_tier_reverse_dict = { item : tier for tier, items in classAfamily_tier_dict.items() for item in items }

# [[[ 1st tier ]]]
classAfamily_color_dict = {}
classAfamily_color_dict["Adrenoceptors"] = "goldenrod"
classAfamily_color_dict["Adenosine"]     = "blue"
classAfamily_color_dict["Opsins"]        = "red"

ps_dict = {"tier 1" : 1.0, "tier 2" : 0.8, "tier 3" : 0.5, "tier 4": 0.3}
## state_dict = {
##                   "Inactive" : {"col" : 11, "shape" : 7, "lw" : 1, },
##                   "Intermediate" : {"col" : 11, "shape" : 7, "lw" : 1, },
##                   "Active"   : {"col" : 11, "shape" : 7, "lw" : 1, }
##              }
is_exist_dict = {}
hl_name = "Class A"
for state_name, attr_dict in state_dict.items():
    color_key = f"{hl_name} {state_name}"
    col = attr_dict["col"]
    shape = attr_dict["shape"]
    lw = attr_dict["lw"]
    line_dict[color_key] = { i : line for i, line in enumerate(lines) if line[col] == state_name }

    # States...
    coloritems[color_key] = select_items(line_dict[color_key], 3)
    colors_dict[color_key] = {**classAfamily_color_dict}
    for family_name, v in colors_dict[color_key].items():
        # Deal with wrongly-input receptor
        if not family_name in coloritems[color_key]:
            print(f"!!! {family_name} is not a valid record for coloring for style '{color_key}'.")
            continue

        # Set title...
        # Title will affect the making of key
        ## t = f"title '{state_name} {family_name}'"

        t = f"title '{family_name}'"

        if not t in is_exist_dict: is_exist_dict[t] = None
        else: t = "notitle"

        ps = ps_dict[classAfamily_tier_reverse_dict[family_name]]
        # Construct plot commands for each entry
        gps_dict[f"{color_key} {family_name}"] = {
            "style" : f"u 1:2 w p pt {shape} ps {ps} lw {lw} lc rgb '{v}' {t}",
            "entry" : coloritems[color_key][family_name]
        }


# [[[ 2nd tier ]]]
classAfamily_color_dict = {}
classAfamily_color_dict.update(cs.color_species(classAfamily_tier_dict["tier 2"], s = 60, v = 90, hexsym = "#", b = 170, e = 340))

## state_dict = {
##                   "Inactive" : {"col" : 11, "shape" : 7, "lw" : 1, },
##                   "Intermediate" : {"col" : 11, "shape" : 7, "lw" : 1, },
##                   "Active"   : {"col" : 11, "shape" : 7, "lw" : 1, }
##              }

is_exist_dict = {}
hl_name = "Class A"
for state_name, attr_dict in state_dict.items():
    color_key = f"{hl_name} {state_name}"
    col = attr_dict["col"]
    shape = attr_dict["shape"]
    lw = attr_dict["lw"]
    line_dict[color_key] = { i : line for i, line in enumerate(lines) if line[col] == state_name }

    # States...
    coloritems[color_key] = select_items(line_dict[color_key], 3)
    colors_dict[color_key] = {**classAfamily_color_dict}
    for family_name, v in colors_dict[color_key].items():
        # Deal with wrongly-input receptor
        if not family_name in coloritems[color_key]:
            print(f"!!! {family_name} is not a valid record for coloring for style '{color_key}'.")
            continue

        # Set title...
        # Title will affect the making of key
        ## t = f"title '{state_name} {family_name}'"

        t = f"title '{family_name}'"

        if not t in is_exist_dict: is_exist_dict[t] = None
        else: t = "notitle"

        ps = ps_dict[classAfamily_tier_reverse_dict[family_name]]
        # Construct plot commands for each entry
        gps_dict[f"{color_key} {family_name}"] = {
            "style" : f"u 1:2 w p pt {shape} ps {ps} lw {lw} lc rgb '{v}' {t}",
            "entry" : coloritems[color_key][family_name]
        }

# [[[ 3rd and 4th tier ]]]
classAfamily_color_dict = {}
classAfamily_color_dict.update(cs.color_species(classAfamily_tier_dict["tier 3"], s = 60, v = 90, hexsym = "#", b = 70, e = 140))
for i in classAfamily_tier_dict["tier 4"]: classAfamily_color_dict[i] = "black"

## state_dict = {
##                   "Inactive" : {"col" : 11, "shape" : 7, "lw" : 1, },
##                   "Intermediate" : {"col" : 11, "shape" : 7, "lw" : 1, },
##                   "Active"   : {"col" : 11, "shape" : 7, "lw" : 1, }
##              }

is_exist_dict = {}
hl_name = "Class A"
for state_name, attr_dict in state_dict.items():
    color_key = f"{hl_name} {state_name}"
    col = attr_dict["col"]
    shape = attr_dict["shape"]
    lw = attr_dict["lw"]
    line_dict[color_key] = { i : line for i, line in enumerate(lines) if line[col] == state_name }

    # States...
    coloritems[color_key] = select_items(line_dict[color_key], 3)
    colors_dict[color_key] = {**classAfamily_color_dict}
    for family_name, v in colors_dict[color_key].items():
        # Deal with wrongly-input receptor
        if not family_name in coloritems[color_key]:
            print(f"!!! {family_name} is not a valid record for coloring for style '{color_key}'.")
            continue

        # Set title...
        # Title will affect the making of key
        ## t = f"title '{state_name} {family_name}'"

        t = f"title '{family_name}'" if classAfamily_tier_reverse_dict[family_name] != "tier 4" else f"title 'Others in family A'"

        if not t in is_exist_dict: is_exist_dict[t] = None
        else: t = "notitle"

        ps = ps_dict[classAfamily_tier_reverse_dict[family_name]]
        # Construct plot commands for each entry
        gps_dict[f"{color_key} {family_name}"] = {
            "style" : f"u 1:2 w p pt {shape} ps {ps} lw {lw} lc rgb '{v}' {t}",
            "entry" : coloritems[color_key][family_name]
        }


cs.color_table(gps_dict)

# [[[ PLOT ]]]
# Create labels...
entry_dict = { f"{v[7]}_{v[10]}" : i for i, v in enumerate(lines) }

# Selectively plot entries...
pdb_list = [ "5NM4_A", "5NLX_A", "6D26_A", "6D27_A",

             "6PS0_A", 

             "6OYA_R", "6LI2_A",

             "4GRV_A", "5WS3_A",

             "6D9H_R",

             "5IUA_A",

             "2YCZ_B",

             "6WQA_A", 

             "2YDO_A",

             "5TE5_A",

             "6PS4_A", "6H7L_A",

             "5OLV_A", "5OLG_A",

             "6PT0_R",
            ]

entry_fil_dict = { k : entry_dict[k] for k in pdb_list if k in entry_dict }
## entry_fil_dict = entry_dict.copy()

cmds = [
           f"set xzeroaxis",
           f"set yzeroaxis",
       ]

top = 41 + 1


def run(rank1, rank2):
    offset = "1.8,0.0"
    ## cmds.append(f"set size ratio -1")
    return plot_coeff(c, rank1, rank2,
                      lbl = gps_dict,
                      label = True, 
                      ## label = False, 
                      label_dict = entry_fil_dict,
                      lbl_fontsize = 6,
                      ## xrange = (-1.5, 1.5),
                      ## yrange = (-1.5, 1.5),

                      ## xrange = (-1.5, 1.5),
                      ## yrange = (-1.0, 2.0),

                      ## xrange = (-0.0, 2.0),
                      ## yrange = (-0.4, -0.1),

                      offset = offset, 
                      rot = 0,
                      height = 2.873,
                      width = 2.873,
                      fontsize = 14,
                      pointsize = 1.0,
                      linewidth = 1.0,
                      fl_path = fl_path,
                      fl_postfix = f'',
                      is_rug = False,
                      index_from_zero = False,
                      cmds = cmds)


if 1:
    if 1:
        # Visualize singular values...
        top = 10
        plot_singular(s, top = top, log = True, 

                      width = 6,
                      height = 5,
                      fontsize = 32,

                      ## width = 4,
                      ## height = 3,
                      ## fontsize = 16,

                      linewidth = 1.0,
                      ticscale  = 2.0,
                      fl_path = fl_path,
                      index_from_zero = False)
    if 0:
        rank1, rank2 = 2, 6
        run(rank1, rank2)

    if 0:
        rank1_adjust, rank2_adjust = rank1_last, rank2_last
        num_job = 10
        start_index = 1
        job_rng = range(1, top)
        ## job_ids = product(job_rng, job_rng)
        job_ids = filter(lambda x: x[0] < x[1], product([rank1_adjust], job_rng))
        with mp.Pool(num_job) as proc:
            proc.starmap( run, job_ids )
        job_ids = filter(lambda x: x[0] < x[1], product([rank2_adjust], job_rng))
        with mp.Pool(num_job) as proc:
            proc.starmap( run, job_ids )
        job_ids = filter(lambda x: x[0] < x[1], product(job_rng, [rank1_adjust]))
        with mp.Pool(num_job) as proc:
            proc.starmap( run, job_ids )
        job_ids = filter(lambda x: x[0] < x[1], product(job_rng, [rank2_adjust]))
        with mp.Pool(num_job) as proc:
            proc.starmap( run, job_ids )

    if 1:
        start_index = 1
        top = 20
        job_rng = range(1, top)
        job_ids = filter(lambda x: x[0] < x[1], permutations(job_rng, 2))
        num_job = 10
        if __name__ == "__main__":
            with mp.Pool(num_job) as proc:
                proc.starmap( run, job_ids )



    if 0:
        offset = "1.8,0.0"
        start_index = 1
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
                           ## xrange = (-2.0, 2.0),
                           ## yrange = (-2.0, 2.0),
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

    if 0:
        offset = "1.8,0.0"
        cmds.append("set size ratio -1")
        rank1, rank2 = 4, 2
        plot_coeff(c, rank1, rank2,
                   lbl = gps_dict,
                   label = True, 
                   ## label_dict = entry_dict,
                   label_dict = entry_fil_dict,
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
