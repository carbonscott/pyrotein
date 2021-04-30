#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyrotein as pr
import numpy as np
import colorsimple as cs
import GnuplotPy3
import os
import tempfile


def plot_dmat(
    dmat,                    # Input data, which is a distance matrix
    fl_dmat,                 # Filename of the exported file
    lbl,                     # Labels used to mark on the diagonal
    lbl_fontsize = 8,        # Fontsize for label
    width         = 6,       # inch
    height        = 7,       # inch
    fontsize      = 14,      # pt
    linewidth     = 1.0,     # pt
    palette       = "",      # Palette definition
    intst_min     = "0",     # Min intensity value
    intst_max     = "*",     # Max intensity value
    vrange        = [],
    showzero      = True,
    showcolorbox  = True,
    NaN           = "NaN",
    mode          = "image", # "image", "sparse", "pm3d"
    cmds_top      = [],      # Customized command for upper panel
    cmds_bottom   = [],      # Customized command for bottom panel
    ):
    assert len(vrange) == 0 or len(vrange) == 2, "vrange has to be an empty or 2-member tuple.  "

    # Partial???
    range_default = ("*", "*")
    if len(vrange) == 2: fl_dmat = f"{fl_dmat}.zoom"

    # Get the mean...
    column_mean_dmat = np.nanmean(dmat, axis = 0, keepdims = False)

    # Draw lbl (optional)...
    cmds_lbl_top = [""]
    cmds_lbl_bottom = [""]
    color_lbl = '#BBBBBB'
    if len(lbl) > 0: 
        for k, (b,e) in lbl.items():
            # Vertical lines (beginning of a region)
            cmd = f"set arrow front from {b-1},graph 0 to {b-1},graph 1 nohead dashtype 2 linewidth {linewidth} linecolor rgb '{color_lbl}'"
            cmds_lbl_bottom.append(cmd)
            cmds_lbl_top.append(cmd)

            # Vertical lines (end of a region)
            cmd = f"set arrow front from {e-1},graph 0 to {e-1},graph 1 nohead dashtype 2 linewidth {linewidth} linecolor rgb '{color_lbl}'"
            cmds_lbl_bottom.append(cmd)
            cmds_lbl_top.append(cmd)

            # Horizontal lines (beginning of a region)
            cmd = f"set arrow front from graph 0,first {b-1} to graph 1,first {b-1} nohead dashtype 2 linewidth {linewidth} linecolor rgb '{color_lbl}'"
            cmds_lbl_bottom.append(cmd)

            # Horizontal lines (end of a region)
            cmd = f"set arrow front from graph 0,first {e-1} to graph 1,first {e-1} nohead dashtype 2 linewidth {linewidth} linecolor rgb '{color_lbl}'"
            cmds_lbl_bottom.append(cmd)

            # Put labels on the diagonal...
            lbl[k] = [ (b + e) // 2, (b + e) // 2 ]

    # [[[ Visualize ]]]
    num_items = len(dmat)
    if intst_max == "*":
        intst_min = np.nanmin(dmat)
        intst_max = np.nanmax(dmat)
    intst_column_mean_min = np.min( [np.nanmin(column_mean_dmat), 0] )
    intst_column_mean_max = np.max( [np.nanmax(column_mean_dmat), 0] )

    # Create tempfile to visualize half matrix...
    fl_temp = tempfile.mktemp(".temp.dat")
    with open(fl_temp,'w') as fh:
        for j in range(num_items):
            for k in range(j if mode != 'image' else num_items):
                if mode == "sparse":
                    if not intst_min < dmat[j, k] < intst_max: continue
                val = NaN if np.isnan(dmat[j, k]) else dmat[j, k]
                fh.write(f"{k} {j} {val}\n")
            fh.write("\n")

    # Begin Gnuplot
    gp = GnuplotPy3.GnuplotPy3()
    gp(f"set terminal postscript eps  size {width}, {height} \\")
    gp(f"                             enhanced color \\")
    gp(f"                             font 'Helvetica,{fontsize}' \\")
    gp(f"                             linewidth {linewidth}")

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
    gp("set bmargin at screen 0.70")
    gp("set lmargin at screen 0.10")
    gp("set rmargin at screen 0.85")
    gp(f"set xrange [-1:{num_items}]")
    gp(f"set yrange [{intst_column_mean_min}:{intst_column_mean_max}]")
    gp("set key top right")
    gp(f"set border linewidth {linewidth}")
    gp("set view map")

    if showzero: gp(f"set arrow front from graph 0, first 0 to graph 1, first 0 nohead dashtype 2 linewidth 1.0 linecolor rgb 'black'")

    for cmd in cmds_lbl_top:
        gp(cmd)

    for cmd in cmds_top:
        gp(cmd)

    gp(f"splot '-' using 1:2:3 with lines linewidth 1.0 linecolor rgb 'black' title 'Column mean'")
    for i,v in enumerate(column_mean_dmat):
        gp(f"{i} {v} 0")
    gp("e")


    # PLOT 2: distance matrix...
    gp(f"unset arrow")
    gp(f"unset key")
    gp(f"unset xrange")
    gp(f"unset yrange")
    gp(f"unset xtics")
    gp(f"unset ytics")
    gp(f"unset logscale")
    gp("unset bmargin")
    gp("unset tmargin")
    gp("unset lmargin")
    gp("unset rmargin")
    gp("set origin 0,0.0")
    gp("set size   1.0,0.70")
    gp("set tmargin at screen 0.7")
    gp("set bmargin at screen 0.05")
    gp("set lmargin at screen 0.10")
    gp("set rmargin at screen 0.85")
    gp(f"set xrange [-1          :{num_items}   ]")
    gp(f"set yrange [{num_items}   :-1          ]")
    gp(f"set border linewidth {linewidth}")

    for k, (x, y) in lbl.items():
        gp(f"set label '{k}' at {x},{y} left rotate by 45 font ', {lbl_fontsize}' front")

    if palette == "":
        gp("set palette defined ( -0.001 'white', 0 'blue', 0.5 'light-grey', 1 'red' )")
    else:
        gp(palette)
    gp(f"set cbrange [{intst_min}:{intst_max}]")
    gp(f"set cbtics font ',{lbl_fontsize}'")
    if not showcolorbox: gp(f"unset colorbox")

    for cmd in cmds_lbl_bottom:
        gp(cmd)

    for cmd in cmds_bottom:
        gp(cmd)

    if mode == 'sparse':
        gp(f"plot '{fl_temp}' using 1:2:3 with points pointtype 6 pointsize 0.5 linewidth 0.5 linecolor palette")
    if mode == 'image':
        gp(f"plot '{fl_temp}' using 1:2:3 with image")
    if mode == 'pm3d':
        gp("set view map")
        gp(f"splot '{fl_temp}' using 1:2:3 with pm3d")

    gp("unset multiplot")
    gp("exit")

    return None




def plot_singular(s, top = 3, fl_export = "singular", 
                              width = 5.65, 
                              height = 5.65,
                              fontsize = 16,
                              linewidth = 1,
                              pointsize = 2,
                              ticscale  = 1.0,
                              fl_path   = '.',
                              log = False, index_from_zero = True):
    ''' Plot singular values.
    '''
    index_base = 0 if index_from_zero else 1

    gp = GnuplotPy3.GnuplotPy3()

    gp(f"set terminal postscript eps  size {width}, {height} \\")
    gp(f"                             enhanced color \\")
    gp(f"                             font 'Helvetica,{fontsize}' \\")
    gp(f"                             linewidth {linewidth}")

    # Set ticscale...
    gp(f"set tic scale {ticscale}")

    # Declare the filename to export...
    fl_export = os.path.join(fl_path, fl_export)
    gp(f"set output '{fl_export}.eps'")
    gp("unset key")

    gp("set xrange [-5:*]")
    if log: gp("set log y")
    gp("set xlabel 'Rank of singular values'")
    gp("set ylabel 'Singular values'")
    gp("plot \\")
    for c in ("#0AFF00", "#000000"):
        gp(f"'-' using 1:2 with points pointsize {pointsize} pointtype 7 linewidth 1 linecolor rgb '{c}' notitle, \\")
    gp("'-' using 1:2 with lines linewidth 1 linecolor rgb '#000000' notitle, \\")
    gp("")

    for i, v in enumerate(s[:top]): 
        gp(f"{i + index_base} {v}")
    gp("e")
    for i, v in enumerate(s[top:]): 
        gp(f"{i+top + index_base} {v}")
    gp("e")
    for i, v in enumerate(s): 
        gp(f"{i + index_base} {v}")
    gp("e")
    gp("exit")

    return None




def plot_left_singular(u, rank, length_mat, 
                                guidelines      = {}, 
                                width           = 6,
                                height          = 7,
                                linewidth       = 1.0,
                                fontsize        = 14,
                                lbl_fontsize    = 10,
                                palette         = None,
                                vrange          = [],
                                frac            = 0.1, 
                                binning         = 4, 
                                fl_path         = '.', 
                                intst_min       = None,
                                intst_max       = None,
                                cmds_top        = [""],
                                cmds_bottom     = [""],
                                fl_postfix      = '',
                                mode            = 'image',
                                index_from_zero = True):
    ''' Plot left singular value as a lower triangular distance matrix.
    '''
    # Comply with the convention (1-based index)
    rank_in_data = rank if index_from_zero else rank - 1

    # Convert `u` at position `rank` into a lower triangular matrix...
    dmat = pr.utils.array2tril(u[:, rank_in_data], length_mat, offset = -1)

    # Restore dmat to a full matrix for visualization...
    dmat_full = dmat + dmat.T

    # Define a color palette...
    # Colorscheme is inspired by [this paper](https://academic.oup.com/nar/article/44/15/7457/2457750)
    if palette is None:
        palette = "set palette defined \
               (-10 '#F6FF9E', -10 '#800000', -5 'red', -1 'white', 0 'seagreen', \
                  1 'white'  , 5 'blue', 10 'navy')"

    # Filename to export...
    fl_export = os.path.join(fl_path, f"u{rank:02d}" + fl_postfix)

    # Bin image???
    dmat_bin = dmat_full
    if binning != 1: dmat_bin = pr.utils.bin_image(dmat, binning = binning)

    # Find the full range...
    bound = np.max(np.abs([np.min(dmat_bin), np.max(dmat_bin)]))
    intst_min = -bound * frac if intst_min == None else intst_min
    intst_max =  bound * frac if intst_max == None else intst_max

    # Visualization...
    plot_dmat(dmat_bin, 
              fl_export, 
              lbl           = guidelines,
              lbl_fontsize =  lbl_fontsize,
              intst_min     = intst_min,
              intst_max     = intst_max,
              width         = width,      # inch
              height        = height,     # inch
              fontsize      = fontsize,
              linewidth     = linewidth,
              palette       = palette, 
              vrange        = vrange,
              mode          = mode,)

    return None




def plot_coeff(c, rank1, rank2, plot_dict = {},
                                label = True,
                                labeltext = [],
                                xrange = ("*", "*"),
                                yrange = ("*", "*"),
                                offset = '2.0,0.0',
                                rot = 0,
                                height = 6,
                                width = 6,
                                fontsize = 16,
                                lbl_fontsize = 4,
                                linewidth = 1.0,
                                pointsize = 1.0,
                                fl_path = '.', 
                                fl_postfix = '',
                                index_from_zero = True,
                                cmds = []):
    ''' Scatter plot of examples from 2 dimensions specified by rank1 and rank2.
    '''
    # Comply with the convention (1-based index)
    rank1_in_data = rank1 if index_from_zero else rank1 - 1
    rank2_in_data = rank2 if index_from_zero else rank2 - 1
    assert rank1_in_data >= 0, "Wrong value for rank1 (base-0 or base-1, is the input index correct?)"
    assert rank2_in_data >= 0, "Wrong value for rank2 (base-0 or base-1, is the input index correct?)"

    if len(plot_dict): 
        gp = GnuplotPy3.GnuplotPy3()
        gp(f"set terminal postscript eps  size {width}, {height} \\")
        gp( "                             enhanced color \\")
        gp(f"                             font 'Helvetica,{fontsize}' \\")
        gp(f"                             linewidth {linewidth}")

        # Declare the filename to export...
        fl_name = f"coeff_{rank1:02d}vs{rank2:02d}" + fl_postfix
        fl_out = os.path.join(fl_path, fl_name)

        ## # Zoom???
        ## range_default = ("*", "*")
        ## if xrange != range_default  or \
        ##    yrange != range_default: fl_out = f"{fl_out}.zoom"

        ## # Label???
        ## if not label: fl_out = f"{fl_out}.nolabel"

        # Decide the final filename
        gp(f"set output '{fl_out}.eps'")
        gp("unset key")

        gp(f"set xrange [{xrange[0]}:{xrange[1]}]")
        gp(f"set yrange [{yrange[0]}:{yrange[1]}]")

        gp(f"set xlabel 'c_{{{rank1:d}}} (\305)'")
        gp(f"set ylabel 'c_{{{rank2:d}}} (\305)'")
        gp("set size 1.0,1.0")
        ## gp("set size ratio -1")

        for cmd in cmds:
            gp(cmd)

        # Plot style...
        gp("plot \\")
        for k, v in plot_dict.items(): gp(" '-' " + v["style"] + ", \\")

        # Label each dot
        if label: gp(f" '-' using 1:2:3:4 with labels rotate variable offset char {offset} font ',{lbl_fontsize}', \\")

        # End plot style
        gp("")

        # Plot entry...
        for k, v in plot_dict.items():
            for i in v["entry"]:
                gp(f"{c[rank1_in_data, i]} {c[rank2_in_data,i]}")
            gp("e")

        if label:
            # Label each dot that is colored only (even thought it's selected from metadata)
            for i, point_label in enumerate(labeltext):
                gp(f"{c[rank1_in_data, i]} {c[rank2_in_data,i]} {point_label} {rot}")  
            gp("e")

        gp("exit")
    else:
        print("Nothing to plot!!!")

    return None




def plot_blankcoeff(rank1, width, height, fl_path, fl_postfix, index_from_zero = True):
    # Comply with the convention (1-based index)
    rank1_in_data = rank1 if index_from_zero else rank1 - 1
    assert rank1_in_data >= 0, "Wrong value for rank1 (base-0 or base-1, is the input index correct?)"

    gp = GnuplotPy3.GnuplotPy3()
    gp(f"set terminal postscript eps  size {width}, {height} \\")
    gp( "                             enhanced color \\")
    gp(f"                             font 'Helvetica,14' \\")
    gp(f"                             linewidth 2")

    # Declare the filename to export...
    fl_name = f"coeff_{rank1:02d}vs{rank1:02d}" + fl_postfix
    fl_out = os.path.join(fl_path, fl_name)

    # Decide the final filename
    gp(f"set output '{fl_out}.eps'")
    gp("unset key")
    gp("set xrange [1:2]")
    gp("set yrange [1:2]")
    gp("unset border")
    gp("unset xtics")
    gp("unset ytics")

    gp("plot \\")
    gp(f"'-' using 1:2 with points pointtype 7 notitle,\\")
    gp("")

    gp(f"0 0")
    gp( "e")
    gp("exit")

    return None




def select_items(lines, col, offset = 0):
    ''' Select the value from column (col).  
    '''
    citems = {}
    for i, v in enumerate(lines):
        val = v[col]
        if not val in citems: citems[val] = [i + offset]
        else: citems[val].append(i + offset)
    return citems




def showHistogram(data, bin_cap, rng, title, cmds = []):
    # Find histogram...
    data_val, data_rng = pr.utils.freq_density(data, bin_cap = bin_cap)

    gp = GnuplotPy3.GnuplotPy3()
    gp("set terminal postscript eps  size 3.5, 2.62 \\")
    gp("                             enhanced color \\")
    gp("                             font 'Helvetica,10' \\")
    gp("                             linewidth 1")
    gp(f"set output '{title}.eps'")
    gp("unset key")

    for cmd in cmds:
        gp(cmd)

    gp("plot '-' using 1:2 with lines linewidth 1 linecolor rgb 'black'")

    for i in range(len(data_val)): 
        if data_rng[i] < rng[0]: continue
        if data_rng[i+1] > rng[1]: continue
        gp(f"{data_rng[i]} {data_val[i]}")  
        gp(f"{data_rng[i+1]} {data_val[i]}")  
    gp("e")

    gp("exit")

    return None
