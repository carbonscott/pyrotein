#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pyrotein as pr
from loaddata import load_xlsx, label_TMs
from display import plot_dmat


# load the data...
rmsd_dmat = np.load("rmsd_dmat.npy")

rmsd_dmat_bin = rmsd_dmat
## rmsd_dmat_bin = pr.utils.bin_image(rmsd_dmat, binning = 1, nan_replace = -1)
fl_rmsd_dmat  = "rmsd"

# Define a colorscheme...
# Colorscheme is inspired by from this publication (DOI: 10.1093/nar/gkw555) from Zhong Ren
pal = "set palette defined ( 0 'white', 0 'seagreen', 0.1 'white', 0.5 'blue', 1 'navy' )"

# Create labels...
labels_TM = label_TMs()
for k, v in labels_TM.items(): labels_TM[k] = [ i * 4 for i in v ]

# Draw guidelines (optional)...
lbls = {}
cmds_guideline_top = [""]
cmds_guideline_bottom = [""]
color_guideline = '#BBBBBB'
if len(labels_TM) > 0: 
    for k, (b,e) in labels_TM.items():
        # Vertical lines (beginning of a region)
        cmd = f"set arrow front from {b-1},graph 0 to {b-1},graph 1 nohead dashtype 2 linewidth 0.5 linecolor rgb '{color_guideline}'"
        cmds_guideline_bottom.append(cmd)
        cmds_guideline_top.append(cmd)

        # Vertical lines (end of a region)
        cmd = f"set arrow front from {e-1},graph 0 to {e-1},graph 1 nohead dashtype 2 linewidth 0.5 linecolor rgb '{color_guideline}'"
        cmds_guideline_bottom.append(cmd)
        cmds_guideline_top.append(cmd)

        # Horizontal lines (beginning of a region)
        cmd = f"set arrow front from graph 0,first {b-1} to graph 1,first {b-1} nohead dashtype 2 linewidth 0.5 linecolor rgb '{color_guideline}'"
        cmds_guideline_bottom.append(cmd)

        # Horizontal lines (end of a region)
        cmd = f"set arrow front from graph 0,first {e-1} to graph 1,first {e-1} nohead dashtype 2 linewidth 0.5 linecolor rgb '{color_guideline}'"
        cmds_guideline_bottom.append(cmd)

        # Put labels on the diagonal...
        lbls[k] = [ (b + e) // 2, (b + e) // 2 ]

# Draw the rigid framework...
fl_fwk = 'fwk.dat'
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


plot_dmat( rmsd_dmat_bin, 
           fl_rmsd_dmat, 
           lbl = lbls,
           width = 10,
           height = 12,
           fontsize = 29,
           lbl_fontsize = 29,
           linewidth = 2.0,
           palette = pal,
           intst_max = "*",
           upper = 0.0,
           smooth = False,
           cmds_top      = cmds_guideline_top,
           cmds_bottom   = cmds_guideline_bottom,
         )


