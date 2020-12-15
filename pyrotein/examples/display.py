#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import GnuplotPy3


def plot_dmat(
    dmat,                  # Input data, which is a distance matrix
    fl_dmat,               # Filename of the exported file
    lbl,                   # Labels used to mark on the diagonal
    lbl_font_size = 8,     # Fontsize for label
    palette       = "",    # Palette definition
    intst_max     = "*",   # Max intensity value
    smooth        = False, # Choices of styles: smooth vs pixelated
    upper         = 0.5,   # The value to facilitate the plot of lower triangle matrix
    cmds_top      = [],    # Customized command for upper panel
    cmds_bottom   = [],    # Customized command for bottom panel
    ):
    # Get the mean...
    column_mean_dmat = np.nanmean(dmat, axis = 0, keepdims = False)

    # [[[ Visualize ]]]
    num_items = len(dmat)
    smooth = False
    if intst_max == "*":
        intst_min = np.nanmin(dmat)
        intst_max = np.nanmax(dmat)
    intst_column_mean_min = np.nanmin(column_mean_dmat)
    intst_column_mean_max = np.nanmax(column_mean_dmat)

    gp = GnuplotPy3.GnuplotPy3()
    gp("set terminal postscript eps  size 6, 7 \\")
    gp("                             enhanced color \\")
    gp("                             font 'Helvetica,14' \\")
    gp("                             linewidth 2")
    gp(f"set output '{fl_dmat}'")
    gp("unset key")

    # Declare a multiplot...
    gp("set origin 0,0")
    gp("set size 1,1")
    gp("unset bmargin")
    gp("unset tmargin")
    gp("unset lmargin")
    gp("unset rmargin")
    gp("set multiplot title ''")


    # PLOT 1: mean dmat...
    gp(f"unset xrange")
    gp(f"unset yrange")
    gp("unset xtics")
    gp("unset ytics")
    gp(f"unset logscale")
    gp("set origin 0,0.70")
    gp("set size   1,0.15")
    gp("set tmargin 0")
    gp("set bmargin 0")
    gp("set lmargin at screen 0.10")
    gp("set rmargin at screen 0.85")
    gp("set border linewidth 0.25")
    gp(f"set xrange [-1:{num_items}]")
    gp(f"set yrange [{intst_column_mean_min}:{intst_column_mean_max}]")
    gp("set key top right")

    for cmd in cmds_top:
        gp(cmd)

    gp(f"plot '-' using 1:2 with linespoints pointtype 6 pointsize 0.4 linewidth 0.25 linecolor rgb 'black' title 'Column mean'")
    for i,v in enumerate(column_mean_dmat):
        gp(f"{i} {v}")
    gp("e")


    # PLOT 2: distance matrix...
    gp(f"unset xrange")
    gp(f"unset yrange")
    gp(f"unset xtics")
    gp(f"unset ytics")
    gp(f"unset logscale")
    gp("set origin 0,0.0")
    gp("set size   1.0,0.70")
    ## gp("set size ratio -1")
    gp("set tmargin 0")
    gp("set bmargin at screen 0.05")
    gp("set lmargin at screen 0.10")
    gp("set rmargin at screen 0.85")
    gp(f"set xrange [-1          :{num_items}   ]")
    gp(f"set yrange [{num_items}   :-1          ]")
    gp("set lmargin at screen 0.10")
    gp("set rmargin at screen 0.85")
    gp("set border linewidth 0.25")

    for j in range(0,num_items):
        gp(f"set label '{lbl[j]}' at {j+0},{j-0} left rotate by 45 font ', {lbl_font_size}' front")

    if palette == "":
        gp("set palette defined ( -0.001 'white', 0 'blue', 0.5 'light-grey', 1 'red' )")
    else:
        gp(palette)
    gp(f"set cbrange [0:{intst_max}]")

    for cmd in cmds_bottom:
        gp(cmd)

    if smooth:
        gp("set pm3d map")
        gp("set pm3d interpolate 0,0")
        gp("set lmargin at screen 0.01")
        gp("set rmargin at screen 0.99")
        gp("set bmargin at screen 0.05")
        gp("set tmargin at screen 0.95")
        gp("splot '-' using 1:2:3")
    else: 
        gp("plot '-' using 1:2:3 with image")
    for j in range(num_items):
        for k in range(num_items):
            if j > k: gp(f"{k} {j} {dmat[j, k]}")
            else: gp(f"{k} {j} {upper}")
        gp(" ")
    gp("e")
    gp("exit")
