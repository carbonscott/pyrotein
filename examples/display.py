#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyrotein as pr
import numpy as np
import colorsimple as cs
import GnuplotPy3
import os
import tempfile


def plot_dmat(
    dmat,                          # Input data, which is a distance matrix
    fl_dmat,                       # Filename of the exported file
    lbl             = {},          # Labels used to mark on the diagonal
    lbl_fontsize    = 8,           # Fontsize for label
    lbl_linewidth   = 1.0,         # pt
    diaglbl         = {},          # diagonal label (usually for showing index)
    diaglblfontsize = 5,
    width           = 6,           # inch
    height          = 7,           # inch
    fontsize        = 14,          # pt
    linewidth       = 1.0,         # pt
    curve_linewidth = 1.0,         # pt
    palette         = "",          # Palette definition
    intst_min       = "0",         # Min intensity value
    intst_max       = "*",         # Max intensity value
    vrange          = [],
    showzero        = True,
    showcolorbox    = True,
    NaN             = "NaN",
    temp            = True,
    mode            = "image",     # "image", "sparse", "pm3d"
    showsparselabel = False,
    box_range       = [],          # If box range is empty, then the box covers the whole area of u
    default_intst_rng = False,
    cmds_top        = [],          # Customized command for upper panel
    cmds_bottom     = [],          # Customized command for bottom panel
    ):
    assert len(vrange) == 0 or len(vrange) == 2, "vrange has to be an empty or 2-member tuple.  "

    # Partial???
    range_default = ("*", "*")
    if len(vrange) == 2: fl_dmat = f"{fl_dmat}.zoom"

    title = f"Column mean"
    cmd_box = [""]
    if len(box_range):
        # Right-inclusion is considered [debating if this should be done]
        b, e = box_range
        # Horizontal lines (beginning of a region)
        cmd = f"set arrow front from graph 0,first {b} to graph 1,first {b} nohead dashtype 2 linewidth {lbl_linewidth} linecolor rgb 'black'"
        cmd_box.append(cmd)

        # Horizontal lines (end of a region)
        cmd = f"set arrow front from graph 0,first {e} to graph 1,first {e} nohead dashtype 2 linewidth {lbl_linewidth} linecolor rgb 'black'"
        cmd_box.append(cmd)

        title += " (within reference range)"
    else:
        b, e = 0, len(dmat)

    # Get the mean...
    column_mean_dmat = np.nanmean(dmat[b:e, :], axis = 0, keepdims = False)

    # Draw lbl (optional)...
    cmds_lbl_top = [""]
    cmds_lbl_bottom = [""]
    color_lbl = '#BBBBBB'
    _lbl = {}
    if len(lbl) > 0: 
        for k, (b,e) in lbl.items():
            # Vertical lines (beginning of a region)
            cmd = f"set arrow front from {b},graph 0 to {b},graph 1 nohead dashtype 2 linewidth {lbl_linewidth} linecolor rgb '{color_lbl}'"
            cmds_lbl_bottom.append(cmd)
            cmds_lbl_top.append(cmd)

            # Vertical lines (end of a region)
            cmd = f"set arrow front from {e},graph 0 to {e},graph 1 nohead dashtype 2 linewidth {lbl_linewidth} linecolor rgb '{color_lbl}'"
            cmds_lbl_bottom.append(cmd)
            cmds_lbl_top.append(cmd)

            # Horizontal lines (beginning of a region)
            cmd = f"set arrow front from graph 0,first {b} to graph 1,first {b} nohead dashtype 2 linewidth {lbl_linewidth} linecolor rgb '{color_lbl}'"
            cmds_lbl_bottom.append(cmd)

            # Horizontal lines (end of a region)
            cmd = f"set arrow front from graph 0,first {e} to graph 1,first {e} nohead dashtype 2 linewidth {lbl_linewidth} linecolor rgb '{color_lbl}'"
            cmds_lbl_bottom.append(cmd)

            # Put labels on the diagonal...
            _lbl[k] = [ (b + e) // 2, (b + e) // 2  ]

    # [[[ Visualize ]]]
    num_items = len(dmat)
    if intst_max == "*":
        intst_min = np.nanmin(dmat)
        intst_max = np.nanmax(dmat)
    intst_column_mean_min = np.min( [np.nanmin(column_mean_dmat), 0] )
    intst_column_mean_max = np.max( [np.nanmax(column_mean_dmat), 0] )

    # Create tempfile to visualize half matrix...
    fl_temp = f"{fl_dmat}.dat"
    if temp: fl_temp = tempfile.mktemp(".temp.dat")
    with open(fl_temp,'w') as fh:
        for j in range(num_items):
            for k in range(j if mode != 'image' else num_items):
                if mode == "sparse":
                    if np.isnan(dmat[j, k]): continue
                    ## if not intst_min < dmat[j, k] < intst_max: continue
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
    gp("set lmargin at screen 0.05")
    gp("set rmargin at screen 0.80")
    gp(f"set xrange [-1:{num_items}]")
    gp(f"set yrange [{intst_column_mean_min}:{1.1 * intst_column_mean_max}]")
    if not default_intst_rng: gp(f"set yrange [{intst_min}:{2.0 * intst_max}]")
    gp("set key top right")
    gp(f"set border linewidth {linewidth}")
    gp("set view map")

    if showzero: gp(f"set arrow front from graph 0, first 0 to graph 1, first 0 nohead dashtype 2 linewidth 1.0 linecolor rgb 'black'")

    for cmd in cmds_lbl_top: gp(cmd)

    for cmd in cmds_top: gp(cmd)

    if mode == "pm3d":
        gp(f"splot '-' using 1:2:3 with lines linewidth {curve_linewidth} linecolor rgb 'black' title '{title}'")
        for i,v in enumerate(column_mean_dmat):
            gp(f"{i} {v} 0")
        gp("e")
    else:
        gp(f"plot '-' using 1:2 with lines linewidth {curve_linewidth} linecolor rgb 'black' title '{title}'")
        for i,v in enumerate(column_mean_dmat):
            gp(f"{i} {v}")
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
    gp("set lmargin at screen 0.05")
    gp("set rmargin at screen 0.80")
    gp(f"set xrange [-1          :{num_items}   ]")
    gp(f"set yrange [{num_items}   :-1          ]")
    gp(f"set border linewidth {linewidth}")

    for k, (x, y) in _lbl.items():
        gp(f"set label '{k}' at {x},{y} left rotate by 45 font ', {lbl_fontsize}' front")

    for k, (x, y) in diaglbl.items():
        gp(f"set label '{k}' at {x},{y} left rotate by 45 font ', {diaglblfontsize}' front")

    if palette == "":
        gp("set palette defined ( -0.001 'white', 0 'blue', 0.5 'light-grey', 1 'red' )")
    else:
        gp(palette)
    gp(f"set cbrange [{intst_min}:{intst_max}]")
    gp(f"set cbtics font ',{lbl_fontsize}'")
    if not showcolorbox: gp(f"unset colorbox")

    for cmd in cmds_lbl_bottom: gp(cmd)

    for cmd in cmds_bottom: gp(cmd)

    for cmd in cmd_box: gp(cmd)

    if mode == 'sparse':
        gp("plot \\")
        gp(f"'{fl_temp}' using 1:2:3 with points pointtype 6 pointsize 0.5 linewidth 0.0 linecolor palette, \\")
        if showsparselabel:
            gp(f"'{fl_temp}' using 1:2:(sprintf('%d,%d', int($1), int($2))) with labels offset 0.5,.3 rotate by 45 font ',3', \\")
        gp("")
    if mode == 'image':
        gp(f"plot '{fl_temp}' using 1:2:3 with image")
    if mode == 'pm3d':
        gp("set view map")
        gp(f"splot '{fl_temp}' using 1:2:3 with pm3d")

    gp("unset multiplot")
    gp("exit")

    return None




def calc_rigid_framework(rmsd_dmat, seqi_fwk_list, num_seq, len_res, min_size = 5, min_mean = 0.1, epsilon = 0.00001):
    ''' The idea of rigid framework is taken from DOI: 10.1371/journal.pone.0077363 .

        Basically, every residue should be assigned to either "rigid" or "not
        rigid".  The change of mean value of the submatrix upon the selection or
        deselection of a residue would determine the assignment.  If the change of
        mean is larger than a threshold, the choice (either select or not) should 
        not be made.  Check residues in both frameworks.  

        idx in rigid fwk, should it be deselected?
        idx not in rigid fwk, should it be selected?
    '''
    atom_list_orig = [ i for b, e in seqi_fwk_list for i in range(b * len_res, e * len_res) ]

    # Calculate the mean of submatrix upon the choice of rigid fwk from input...
    def calc_submean(rmsd_dmat, atom_list_orig):
        rmsd_dmat_sub_aux_orig  = np.take(rmsd_dmat, atom_list_orig, axis = 1)
        rmsd_dmat_sub_orig      = np.take(rmsd_dmat_sub_aux_orig, atom_list_orig, axis = 0)
        mean_rmsd_dmat_sub_orig = np.nanmean(rmsd_dmat_sub_orig, keepdims = False)
        return mean_rmsd_dmat_sub_orig

    mean_rmsd_dmat_sub_orig = calc_submean(rmsd_dmat, atom_list_orig)
    print(f"RMSD (Init): {mean_rmsd_dmat_sub_orig}")

    # Fetch index of framework-wise np.nan...
    mean_rmsd_dmat_orig = np.nanmean(rmsd_dmat, axis = 0, keepdims = False)
    mean_rmsd_dmat_orig[np.isnan(mean_rmsd_dmat_orig)] = 0.0
    nan_list = np.argwhere(mean_rmsd_dmat_orig < min_mean).reshape(-1)

    # Inspect every residue (col) represented by seqi...
    # Easier to operate on atom list, this is a adhoc code anyway
    # - idx in rigid fwk, should it be deselected?
    # - idx not in rigid fwk, should it be selected?
    rm_list  = []
    add_list = []
    for seqi in range(num_seq):
        # Infer atomi from seqi by a scale of len_res...
        atomi = seqi * len_res

        # If this residue is rigid???
        if atomi in atom_list_orig:
            # Remove nan column...
            if atomi in nan_list:
                rm_list.append(atomi)
                continue

            atom_list_aux = atom_list_orig.copy()

            # Remove the associated atoms...
            for i in range(len_res): atom_list_aux.remove(atomi + i)

            # Calcualte the new mean rmsd...
            mean_rmsd_dmat_sub_aux = calc_submean(rmsd_dmat, atom_list_aux)

            # Consider NOT to keep it if it's larger than 1% of decrease???
            if mean_rmsd_dmat_sub_aux > (1 - epsilon) * mean_rmsd_dmat_sub_orig: continue

            rm_list.append(atomi)

        # Or this residue is not rigid???
        else:
            # Don't consider nan...
            if atomi in nan_list: continue

            atom_list_aux = atom_list_orig.copy()

            # Remove the associated atoms...
            for i in range(len_res): atom_list_aux.append(atomi + i)

            # Calcualte the new mean rmsd...
            mean_rmsd_dmat_sub_aux = calc_submean(rmsd_dmat, atom_list_aux)

            # Consider NOT to keep it if it's larger than 1% of increase???
            if mean_rmsd_dmat_sub_aux > (1 + epsilon) * mean_rmsd_dmat_sub_orig: continue

            add_list.append(atomi)

    # Remove atoms that was considered rigid...
    for atomi in rm_list:
        # Remove the associated atoms...
        for i in range(len_res): atom_list_orig.remove(atomi + i)

    # Add atoms that was considered not rigid...
    for atomi in add_list:
        # Remove the associated atoms...
        for i in range(len_res): atom_list_orig.append(atomi + i)

    atom_list_orig.sort()
    fwk_list = pr.utils.group_consecutive_integer(atom_list_orig)
    fwk_list = [ i for i in fwk_list if i[-1] - i[0] + 1 >= min_size * len_res ]
    atom_list_orig = [ i for fwk in fwk_list for i in fwk ]

    mean_rmsd_dmat_sub_orig = calc_submean(rmsd_dmat, atom_list_orig)
    print(f"RMSD (Final): {mean_rmsd_dmat_sub_orig}")

    ## return [ i for i in fwk_list if i[-1] - i[0] + 1 >= min_size * len_res ]
    return fwk_list




def plot_rmsd_dmat(
    dmat,                    # Input data, which is a distance matrix
    fl_dmat,                 # Filename of the exported file
    lbl = {},                # Labels used to mark on the diagonal
    lbl_fontsize = 8,        # Fontsize for label
    diaglbl       = {},      # diagonal label (usually for showing index)
    diaglblfontsize = 5,
    fwk_list = [],
    fwk_linewidth = 2,
    fwk_curve_color = "black",
    fwk_box_color = "black",
    width         = 6,       # inch
    height        = 7,       # inch
    fontsize      = 14,      # pt
    linewidth     = 1.0,     # pt
    lbl_linewidth = 2.0,     # pt
    curve_linewidth = 2.0,     # pt
    curve_color   = "gray",
    palette       = "",      # Palette definition
    intst_min     = "0",     # Min intensity value
    intst_max     = "*",     # Max intensity value
    vrange        = [],
    showzero      = True,
    showcolorbox  = True,
    NaN           = "NaN",
    temp          = True,
    mode          = "image", # "image", "sparse", "pm3d"
    showsparselabel = False,
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
    _lbl = {}
    if len(lbl) > 0: 
        for k, (b,e) in lbl.items():
            # Vertical lines (beginning of a region)
            cmd = f"set arrow front from {b},graph 0 to {b},graph 1 nohead dashtype 2 linewidth {lbl_linewidth} linecolor rgb '{color_lbl}'"
            cmds_lbl_bottom.append(cmd)
            cmds_lbl_top.append(cmd)

            # Vertical lines (end of a region)
            cmd = f"set arrow front from {e},graph 0 to {e},graph 1 nohead dashtype 2 linewidth {lbl_linewidth} linecolor rgb '{color_lbl}'"
            cmds_lbl_bottom.append(cmd)
            cmds_lbl_top.append(cmd)

            # Horizontal lines (beginning of a region)
            cmd = f"set arrow front from graph 0,first {b} to graph 1,first {b} nohead dashtype 2 linewidth {lbl_linewidth} linecolor rgb '{color_lbl}'"
            cmds_lbl_bottom.append(cmd)

            # Horizontal lines (end of a region)
            cmd = f"set arrow front from graph 0,first {e} to graph 1,first {e} nohead dashtype 2 linewidth {lbl_linewidth} linecolor rgb '{color_lbl}'"
            cmds_lbl_bottom.append(cmd)

            # Put labels on the diagonal...
            _lbl[k] = [ (b + e) // 2, (b + e) // 2  ]

    # [[[ Visualize ]]]
    num_items = len(dmat)
    if intst_max == "*":
        intst_min = np.nanmin(dmat)
        intst_max = np.nanmax(dmat)
    intst_column_mean_min = np.min( [np.nanmin(column_mean_dmat), 0] )
    intst_column_mean_max = np.max( [np.nanmax(column_mean_dmat), 0] )

    # Create tempfile to visualize half matrix...
    fl_temp = f"{fl_dmat}.dat"
    if temp: fl_temp = tempfile.mktemp(".temp.dat")
    with open(fl_temp,'w') as fh:
        for j in range(num_items):
            for k in range(j if mode != 'image' else num_items):
                if mode == "sparse":
                    if np.isnan(dmat[j, k]): continue
                    ## if not intst_min < dmat[j, k] < intst_max: continue
                val = NaN if np.isnan(dmat[j, k]) else dmat[j, k]
                fh.write(f"{k} {j} {val}\n")
            fh.write("\n")


    # [[[ RIGID FRAMEWORK ]]]
    for i in range(len(fwk_list)):
        cmd = []
        b1, e1 = [ fwk_list[i][0], fwk_list[i][-1] ]
        print(f"{b1//4} {e1//4}")
        cmd.append(f"set object rectangle front from {b1},{e1} to {e1},{b1} fs empty border linecolor rgb 'black'")
        for j in range(i + 1, len(fwk_list)):
            b2, e2 = [ fwk_list[j][0], fwk_list[j][-1] ]
            cmd.append(f"set object rectangle front from {b1},{b2} to {e1},{e2} fs empty linecolor rgb '{fwk_box_color}'")
            cmd.append(f"set object rectangle front from {b2},{b1} to {e2},{e1} fs empty linecolor rgb '{fwk_box_color}'")
        cmds_bottom.extend(cmd)


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

    if showzero: gp(f"set arrow front from graph 0, first 0 to graph 1, first 0 nohead dashtype 2 linewidth {lbl_linewidth} linecolor rgb 'black'")

    for cmd in cmds_lbl_top:
        gp(cmd)

    for cmd in cmds_top:
        gp(cmd)

    if mode == "pm3d":
        gp(f"splot \\")
        gp(f"'-' using 1:2:3 with lines linewidth {curve_linewidth} linecolor rgb '{curve_color}' title 'Column mean', \\")
        if len(fwk_list): gp(f"'-' using 1:2:3 with lines linewidth {fwk_linewidth} linecolor rgb '{fwk_curve_color}' notitle, \\")
        gp(f"")

        for i,v in enumerate(column_mean_dmat):
            gp(f"{i} {v} 0")
        gp("e")

        for fwk in fwk_list:
            for i in fwk:
                v = column_mean_dmat[i]
                gp(f"{i} {v}")
            gp(" ")
        if len(fwk_list): gp("e")
    else:
        gp(f"plot \\")
        gp(f"'-' using 1:2 with lines linewidth {curve_linewidth} linecolor rgb '{curve_color}' title 'Column mean', \\")
        if len(fwk_list): gp(f"'-' using 1:2 with lines linewidth {fwk_linewidth} linecolor rgb '{fwk_curve_color}' notitle, \\")
        gp("")

        for i,v in enumerate(column_mean_dmat):
            gp(f"{i} {v}")
        gp("e")

        for fwk in fwk_list:
            for i in fwk:
                v = column_mean_dmat[i]
                gp(f"{i} {v}")
            gp(" ")
        if len(fwk_list): gp("e")


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

    for k, (x, y) in _lbl.items():
        gp(f"set label '{k}' at {x},{y} left rotate by 45 font ', {lbl_fontsize}' front")

    for k, (x, y) in diaglbl.items():
        gp(f"set label '{k}' at {x},{y} left rotate by 45 font ', {diaglblfontsize}' front")

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
        gp("plot \\")
        gp(f"'{fl_temp}' using 1:2:3 with points pointtype 6 pointsize 0.5 linewidth 0.0 linecolor palette, \\")
        if showsparselabel:
            gp(f"'{fl_temp}' using 1:2:(sprintf('%d,%d', int($1), int($2))) with labels offset 0.5,.3 rotate by 45 font ',3', \\")
        gp("")
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
                                lbl             = {}, 
                                diaglbl         = {},
                                diaglblfontsize = 5,
                                width           = 6,
                                height          = 7,
                                linewidth     = 1.0,     # pt
                                lbl_linewidth = 2.0,     # pt
                                curve_linewidth = 2.0,     # pt
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
                                temp            = True,
                                showsparselabel = False,
                                default_intst_rng = False,
                                mode            = 'image',
                                box_range       = [],
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
               (-10 '#800000', \
                 -5 'red', \
                 -1 'white', \
                  0 'seagreen', \
                  1 'white'  , \
                  5 'blue', \
                  10 'navy')"

    # Filename to export...
    fl_eps = f"u{rank:02d}"
    if len(box_range): fl_eps += ".box"
    fl_eps += fl_postfix
    fl_export = os.path.join(fl_path, fl_eps)

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
              lbl             = lbl,
              lbl_fontsize    = lbl_fontsize,
              diaglbl         = diaglbl,
              diaglblfontsize = diaglblfontsize,
              intst_min       = intst_min,
              intst_max       = intst_max,
              width           = width,      # inch
              height          = height,     # inch
              fontsize        = fontsize,
              linewidth       = linewidth,
              lbl_linewidth   = lbl_linewidth,
              curve_linewidth = curve_linewidth,
              palette         = palette, 
              vrange          = vrange,
              temp            = temp,
              showsparselabel = showsparselabel,
              default_intst_rng = default_intst_rng,
              box_range       = box_range,
              cmds_top        = cmds_top,
              cmds_bottom     = cmds_bottom,
              mode            = mode,)

    return None




def plot_coeff(c, rank1, rank2, lbl = {},
                                label = True,
                                label_dict = {},
                                join_dict = {},
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
                                precision = '%.1f',
                                fl_path = '.', 
                                fl_postfix = '',
                                is_rug = False,
                                index_from_zero = True,
                                cmds = []):
    ''' Scatter plot of examples from 2 dimensions specified by rank1 and rank2.
    '''
    # Comply with the convention (1-based index)
    rank1_in_data = rank1 if index_from_zero else rank1 - 1
    rank2_in_data = rank2 if index_from_zero else rank2 - 1
    assert rank1_in_data >= 0, "Wrong value for rank1 (base-0 or base-1, is the input index correct?)"
    assert rank2_in_data >= 0, "Wrong value for rank2 (base-0 or base-1, is the input index correct?)"

    if len(lbl): 
        gp = GnuplotPy3.GnuplotPy3()
        gp(f"set terminal postscript eps  size {width}, {height} \\")
        gp( "                             enhanced color \\")
        gp(f"                             font 'Helvetica,{fontsize}' \\")
        gp(f"                             linewidth {linewidth}")
        gp(f"set encoding utf8")

        # Keep the scale on both direction the same...

        # Declare the filename to export...
        fl_name = f"coeff_{rank1:02d}vs{rank2:02d}" + fl_postfix
        fl_out = os.path.join(fl_path, fl_name)

        # Decide the final filename
        gp(f"set output '{fl_out}.eps'")
        gp("unset key")

        gp(f"set xrange [{xrange[0]}:{xrange[1]}]")
        gp(f"set x2range [{xrange[0]}:{xrange[1]}]")    # To facilitate rug plot
        gp(f"set yrange [{yrange[0]}:{yrange[1]}]")
        gp(f"set y2range [{yrange[0]}:{yrange[1]}]")
        gp(f"set format x '{precision}'")
        gp(f"set format y '{precision}'")

        gp(f"set xlabel 'c_{{{rank1:d}}} (\305)'")
        gp(f"set ylabel 'c_{{{rank2:d}}} (\305)'")

        # Rug plot
        if is_rug:
            # Set up tics for rug plot...
            gp("set xtics nomirror")
            gp("set ytics nomirror")
            gp("set x2tics scale 1")
            gp("set y2tics scale 1")
            ## gp("set border lw 0.25")

        for cmd in cmds:
            gp(cmd)

        # Plot style...
        gp("plot \\")
        for k, v in lbl.items(): gp(" '-' " + v["style"] + ", \\")

        # Label each dot
        if label: gp(f" '-' using 1:2:3:4 with labels rotate variable offset char {offset} font ',{lbl_fontsize}', \\")

        # Join dots???
        if bool(join_dict): gp(" '-' using 1:2 with lines linewidth 1.0 linecolor rgb 'gray', \\")

        if is_rug: 
            # Add rug plot...
            for k, v in lbl.items(): 
                gp(" '-' using 1:2:x2tic(''):y2tic('') with points pointsize 0 linewidth 0 notitle, \\")
                gp(" '-' using 1:2:x2tic(''):y2tic('') with points pointsize 0 linewidth 0 notitle, \\")

        # End plot style
        gp("")

        # Plot entry...
        for k, v in lbl.items():
            for i in v["entry"]:
                gp(f"{c[rank1_in_data, i]} {c[rank2_in_data,i]}")
            gp("e")

        if label:
            # Label each dot that is colored only (even thought it's selected from metadata)
            for point_label, i in label_dict.items():
                gp(f"{c[rank1_in_data, i]} {c[rank2_in_data,i]} {point_label} {rot}")  
            gp("e")

        if bool(join_dict):
            for k, v in join_dict.items():
                if len(v) != 2: continue
                for i in v:
                    gp(f"{c[rank1_in_data, i]} {c[rank2_in_data,i]}")  
                gp("")
            gp("e")

        if is_rug:
            # ...w/ rug plot
            for k, v in lbl.items():
                for i in v["entry"]:
                    gp(f"{c[rank1_in_data, i]} {c[rank2_in_data,i]}")
                gp("e")
                for i in v["entry"]:
                    gp(f"{c[rank1_in_data, i]} {c[rank2_in_data,i]}")
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
    gp("set encoding utf8")

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




## def select_items(lines, col, offset = 0):
##     ''' Select the value from column (col).  
##     '''
##     citems = {}
##     for i, v in enumerate(lines):
##         val = v[col]
##         if not val in citems: citems[val] = [i + offset]
##         else: citems[val].append(i + offset)
##     return citems




def select_items(line_dict, col, offset = 0):
    ''' Select the value from column (col).  
    '''
    citems = {}
    for i, line in line_dict.items():
        val = line[col]
        if not val in citems: citems[val] = [i + offset]
        else: citems[val].append(i + offset)
    return citems




def showHistogram(data, bin_cap, rng, title, cmds = []):
    # Find histogram...
    data_val, data_rng = pr.utils.population_density(data, bin_cap = bin_cap)

    gp = GnuplotPy3.GnuplotPy3()
    gp("set terminal postscript eps  size 3.5, 2.62 \\")
    gp("                             enhanced color \\")
    gp("                             font 'Helvetica,10' \\")
    gp("                             linewidth 1")
    gp(f"set output '{title}.eps'")
    gp(f"set encoding utf8")
    gp("unset key")
    gp("set xlabel 'Distance (\305)'")
    gp("set ylabel 'Population density (1/\305)'")
    gp(f"set xrange [{rng[0]}:{rng[1]}]")

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




def plot_u_ave(
    u,                    # Input data, which is a distance matrix
    rank,
    length_mat,
    lbl = {},                # Labels used to mark on the diagonal
    lbl_fontsize = 8,        # Fontsize for label
    xrange = ("*", "*"),
    xstep  = 0,
    yrange = ("*", "*"),
    width         = 5,       # inch
    height        = 1.5,     # inch
    fontsize      = 14,      # pt
    linewidth     = 1.0,     # pt
    intst_min     = "0",     # Min intensity value
    intst_max     = "*",     # Max intensity value
    showzero      = True,
    fl_path       = '.',
    fl_postfix      = '',
    cmds          = [],
    index_from_zero = True,
    ):
    # Comply with the convention (1-based index)
    rank_in_data = rank if index_from_zero else rank - 1

    # Convert `u` at position `rank` into a lower triangular matrix...
    dmat = pr.utils.array2tril(u[:, rank_in_data], length_mat, offset = -1)

    # Restore dmat to a full matrix for visualization...
    dmat_full = dmat + dmat.T

    # Get the mean...
    column_mean_dmat = np.nanmean(dmat_full, axis = 0, keepdims = False)

    # Draw lbl (optional)...
    cmds_lbl_top = [""]
    cmds_lbl_bottom = [""]
    color_lbl = '#BBBBBB'
    if len(lbl) > 0: 
        for k, (b,e) in lbl.items():
            # Vertical lines (beginning of a region)
            cmd = f"set arrow front from {b},graph 0 to {b},graph 1 nohead dashtype 2 linewidth {linewidth} linecolor rgb '{color_lbl}'"
            cmds_lbl_bottom.append(cmd)
            cmds_lbl_top.append(cmd)

            # Vertical lines (end of a region)
            cmd = f"set arrow front from {e},graph 0 to {e},graph 1 nohead dashtype 2 linewidth {linewidth} linecolor rgb '{color_lbl}'"
            cmds_lbl_bottom.append(cmd)
            cmds_lbl_top.append(cmd)

            # Horizontal lines (beginning of a region)
            cmd = f"set arrow front from graph 0,first {b} to graph 1,first {b} nohead dashtype 2 linewidth {linewidth} linecolor rgb '{color_lbl}'"
            cmds_lbl_bottom.append(cmd)

            # Horizontal lines (end of a region)
            cmd = f"set arrow front from graph 0,first {e} to graph 1,first {e} nohead dashtype 2 linewidth {linewidth} linecolor rgb '{color_lbl}'"
            cmds_lbl_bottom.append(cmd)

            # Put labels on the diagonal...
            lbl[k] = (b + e) // 2

    # [[[ Visualize ]]]
    num_items = len(dmat)
    if intst_max == "*":
        intst_min = np.nanmin(dmat)
        intst_max = np.nanmax(dmat)
    intst_column_mean_min = np.min( [np.nanmin(column_mean_dmat), 0] )
    intst_column_mean_max = np.max( [np.nanmax(column_mean_dmat), 0] )

    # Filename to export...
    fl_export = os.path.join(fl_path, f"ave.{rank:02d}{fl_postfix}")

    # Begin Gnuplot
    gp = GnuplotPy3.GnuplotPy3()
    gp(f"set terminal postscript eps  size {width}, {height} \\")
    gp(f"                             enhanced color \\")
    gp(f"                             font 'Helvetica,{fontsize}' \\")
    gp(f"                             linewidth {linewidth}")

    # Declare the filename to export...
    gp(f"set output '{fl_export}.eps'")
    gp("unset key")

    gp(f"set xrange [-1:{num_items}]")
    gp(f"set yrange [{intst_column_mean_min}:{intst_column_mean_max}]")
    gp(f"set border linewidth {linewidth}")

    if xrange != ("*", "*"): gp(f"set xrange [{xrange[0]}:{xrange[1]}]")
    if yrange != ("*", "*"): gp(f"set yrange [{yrange[0]}:{yrange[1]}]")

    for k, x in lbl.items():
        gp(f"set label '{k}' at {x},graph 0.9 left rotate by 90 font ', {lbl_fontsize}' front")

    if showzero: gp(f"set arrow front from graph 0, first 0 to graph 1, first 0 nohead dashtype 2 linewidth 1.0 linecolor rgb 'black'")

    gp(f"set xlabel 'seqi'")
    gp(f"set xtic nomirror")

    if xrange != ("*", "*"):
        b, e = int(xrange[0]), int(xrange[1])
        if not xstep: print(f"!!! xstep must be larger than 0, labeling has failed.")
        else:
            x2l = [ f"'{ i // 4 }_{ i % 4 }' {i}" for i in range(int(b), int(e) + 1, xstep) ]
            cmd_x2l = ', '.join(x2l)
            gp(f"set x2tics ({cmd_x2l}) rotate by 90 right")

    for cmd in cmds_lbl_top:
        gp(cmd)

    for cmd in cmds:
        gp(cmd)

    gp(f"plot '-' using 1:2 with lines linewidth {linewidth} linecolor rgb 'black' title 'Column mean'")
    for i,v in enumerate(column_mean_dmat):
        gp(f"{i} {v}")
    gp("e")

    gp("exit")

    return None
