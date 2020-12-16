#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyrotein as pr
import numpy as np
import GnuplotPy3


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




def plot_singular(s, eps = "singular.eps", log = False):
    ''' Plot singular values.
    '''
    gp = GnuplotPy3.GnuplotPy3()
    gp("set terminal postscript eps size 3.5, 2.62 enhanced color font 'Helvetica,14' linewidth 2")
    gp(f"set output '{eps}'")
    gp("unset key")
    if log: gp("set log y")
    gp("set xlabel 'Rank of singular values'")
    gp("set ylabel 'Singular values'")
    gp("plot '-' using 1:2 with linespoints linetype 2 linewidth 1 linecolor rgb 'black'")

    for i, v in enumerate(s): gp(f"{i} {v}")
    gp("e")
    gp("exit")




def plot_left_singular(u, rank, length_mat, frac = 0.1, binning = 4):
    ''' Plot left singular value as a lower triangular distance matrix.
    '''
    # Convert `u` at position `rank` into a lower triangular matrix...
    dmat = pr.utils.array2tril(u[:, rank], length_mat, offset = -1)

    # Define a color palette...
    # Colorscheme is inspired by [this paper](https://academic.oup.com/nar/article/44/15/7457/2457750)
    pal = "set palette defined \
           (-10.001 'white', -10 '#800000', -5 'red', -1 'white', 0 'seagreen', \
              1 'white'  , 5 'blue', 10 'navy')"

    # Filename to export...
    eps = f"u{rank:02d}.eps"

    # Export eps...
    dmat_bin = pr.utils.bin_image(dmat, binning = binning)
    dmat_bin_check = dmat_bin.copy()

    # Find the full range...
    bound = np.max(np.abs([np.min(dmat_bin), np.max(dmat_bin)]))
    intst_min = -bound * frac
    intst_max =  bound * frac

    # Visualization...
    plot_dmat(dmat_bin, 
              eps, 
              lbl = [""] * len(dmat_bin), 
              intst_min = intst_min,
              intst_max = intst_max,
              palette = pal, 
              upper  = intst_min - 1,
              smooth = True)




def plot_coeff(c, rank1, rank2, labels, offset = '2.0,0.0', rot = 0):
    ''' Scatter plot of examples from 2 dimensions specified by rank1 and rank2.
    '''
    gp = GnuplotPy3.GnuplotPy3()
    gp("set terminal postscript eps size 3.5, 2.62 enhanced color font 'Helvetica,14' linewidth 2")
    gp(f"set output 'coeff_{rank1:02d}vs{rank2:02d}.eps'")
    gp("unset key")
    gp(f"set xlabel 'c_{{{rank1:02d}}}'")
    gp(f"set ylabel 'c_{{{rank2:02d}}}'")

    gp( "plot '-' using 1:2   with point linewidth 1.5 pointtype 6 pointsize 0.5 linecolor rgb '#003f5c', \\")
    gp(f"     '-' using 1:2:3 with labels rotate by {rot} offset char {offset} font ',8', \\")
    gp("")
    for i in range(len(c[0])): 
        gp(f"{c[rank1, i]} {c[rank2,i]}")  
    gp("e")
    for i in range(len(c[0])): 
        gp(f"{c[rank1, i]} {c[rank2,i]} {labels[i]}")  
    gp("e")
    gp("exit")
