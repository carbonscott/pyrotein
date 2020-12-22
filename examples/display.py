#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyrotein as pr
import numpy as np
import colorsimple as cs
import GnuplotPy3
import random
import math


def plot_dmat(
    dmat,                  # Input data, which is a distance matrix
    fl_dmat,               # Filename of the exported file
    lbl,                   # Labels used to mark on the diagonal
    lbl_font_size = 8,     # Fontsize for label
    palette       = "",    # Palette definition
    intst_min     = "0",   # Min intensity value
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
    gp("                             linewidth 1")

    # Declare the filename to export...
    gp(f"set output '{fl_dmat}.eps'")
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

    gp(f"plot '-' using 1:2 with lines linewidth 0.25 linecolor rgb 'black' title 'Column mean'")
    ## gp(f"plot '-' using 1:2 with linespoints pointtype 6 pointsize 0.4 linewidth 0.25 linecolor rgb 'black' title 'Column mean'")
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
    gp(f"set cbrange [{intst_min}:{intst_max}]")

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




def plot_singular(s, top = 3, fl_export = "singular", log = False):
    ''' Plot singular values.
    '''
    gp = GnuplotPy3.GnuplotPy3()

    gp("set terminal postscript eps  size 3.5, 2.62 \\")
    gp("                             enhanced color \\")
    gp("                             font 'Helvetica,14' \\")
    gp("                             linewidth 1.5")

    # Declare the filename to export...
    gp(f"set output '{fl_export}.eps'")
    gp("unset key")

    gp("set xrange [-5:*]")
    if log: gp("set log y")
    gp("set xlabel 'Rank of singular values'")
    gp("set ylabel 'Singular values'")
    gp("plot \\")
    for c in ("#0AFF00", "#000000"):
        gp(f"'-' using 1:2 with points pointtype 7 linewidth 1 linecolor rgb '{c}' notitle, \\")
    gp("'-' using 1:2 with lines linewidth 1 linecolor rgb '#000000' notitle, \\")
    gp("")

    for i, v in enumerate(s[:top]): 
        gp(f"{i} {v}")
    gp("e")
    for i, v in enumerate(s[top:]): 
        gp(f"{i+top} {v}")
    gp("e")
    for i, v in enumerate(s): 
        gp(f"{i} {v}")
    gp("e")
    gp("exit")




def plot_left_singular(u, rank, length_mat, frac = 0.1, binning = 4):
    ''' Plot left singular value as a lower triangular distance matrix.
    '''
    # Convert `u` at position `rank` into a lower triangular matrix...
    dmat = pr.utils.array2tril(u[:, rank], length_mat, offset = -1)

    # Restore dmat to a full matrix for visualization...
    dmat_full = dmat + dmat.T

    # Define a color palette...
    # Colorscheme is inspired by [this paper](https://academic.oup.com/nar/article/44/15/7457/2457750)
    pal = "set palette defined \
           (-10.001 'white', -10 '#800000', -5 'red', -1 'white', 0 'seagreen', \
              1 'white'  , 5 'blue', 10 'navy')"

    # Filename to export...
    fl_export = f"u{rank:02d}"

    # Bin image???
    dmat_bin = dmat_full
    if binning != 1: dmat_bin = pr.utils.bin_image(dmat, binning = binning)

    # Find the full range...
    bound = np.max(np.abs([np.min(dmat_bin), np.max(dmat_bin)]))
    intst_min = -bound * frac
    intst_max =  bound * frac

    # Visualization...
    plot_dmat(dmat_bin, 
              fl_export, 
              lbl = [""] * len(dmat_bin), 
              intst_min = intst_min,
              intst_max = intst_max,
              palette = pal, 
              upper  = intst_min - 1,
              smooth = True)




def plot_coeff(c, rank1, rank2, point_labels,
                                color_items,
                                xrange = ("*", "*"),
                                yrange = ("*", "*"),
                                offset = '2.0,0.0',
                                rot = 0,
                                height = 3,
                                width = 3):
    ''' Scatter plot of examples from 2 dimensions specified by rank1 and rank2.
    '''
    gp = GnuplotPy3.GnuplotPy3()
    gp(f"set terminal postscript eps  size {width}, {height} \\")
    gp( "                             enhanced color \\")
    gp( "                             font 'Helvetica,14' \\")
    gp( "                             linewidth 1.5")

    # Declare the filename to export...
    gp(f"set output 'coeff_{rank1:02d}vs{rank2:02d}.eps'")
    gp("unset key")

    gp(f"set xrange [{xrange[0]}:{xrange[1]}]")
    gp(f"set yrange [{yrange[0]}:{yrange[1]}]")

    gp(f"set xlabel 'c_{{{rank1:02d}}} (\305)'")
    gp(f"set ylabel 'c_{{{rank2:02d}}} (\305)'")
    gp("set size 1.0,1.0")
    ## gp("set size ratio -1")

    gp("plot \\")

    # Generate biolerplate code for Gnuplot
    spe = { i : 0 for i in color_items }.keys()
    ## spe = set(color_items)
    color_dict = cs.color_species(spe)
    for name, color in color_dict.items():
        gp(f"'-' using 1:2   with point linewidth 1.0 pointtype 6 pointsize 0.8 linecolor rgb '{color}' title '{name}', \\")

    # Label each dot
    gp(f"'-' using 1:2:3:4 with labels rotate variable offset char {offset} font ',8', \\")

    # Connecting dots
    gp(f"'-' using 1:2 with lines linewidth 0.5 linecolor rgb '#999999'")
    gp("")

    # The plot statement addressing plot by colors
    for k in color_dict.keys():
        for i in range(len(c[0])): 
            if color_items[i] == k:
                gp(f"{c[rank1, i]} {c[rank2,i]}")  
        gp("e")

    # Label each dot
    ## last_pdb = ""
    ## for i in range(len(c[0])): 
    ##     point_label = point_labels[i]
    ##     curr_pdb = point_label[:4]
    ##     if last_pdb == curr_pdb: point_label = point_label[-1]
    ##     gp(f"{c[rank1, i]} {c[rank2,i]} {point_label}")  
    ##     last_pdb = curr_pdb
    ## gp("e")
    for i in range(len(c[0])): 
        point_label = point_labels[i]
        ## gp(f"{c[rank1, i]} {c[rank2,i]} {point_label} {random.random() * 360}")  
        gp(f"{c[rank1, i]} {c[rank2,i]} {point_label} {rot}")  
    gp("e")

    # Connecting dots
    last_pdb = point_labels[0][:4]
    for i in range(len(c[0])): 
        curr_pdb = point_labels[i][:4]
        if last_pdb != curr_pdb: gp("")
        gp(f"{c[rank1, i]} {c[rank2,i]}")
        last_pdb = curr_pdb
    gp("exit")

    cs.color_table(color_dict)
