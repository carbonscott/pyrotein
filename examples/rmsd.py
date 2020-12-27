#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pyrotein as pr
from loaddata import load_xlsx, label_TMs
from display import plot_dmat


# Specify chains to process...
fl_chain = "chains.comp.xlsx"
lines    = load_xlsx(fl_chain)
drc      = "pdb"

# Define atoms used for distance matrix analysis...
peptide = ["N", "CA", "C", "O"]

# Specify the range of atoms from rhodopsin...
nterm = 1
cterm = 322
len_peptide = (cterm - nterm + 1) * len(peptide)

dmats = np.zeros((len(lines), len_peptide, len_peptide))
for i_fl, line in enumerate(lines):
    # Unpack parameters
    _, pdb, chain, species = line[:4]

    # Read coordinates from a PDB file...
    fl_pdb    = f"{pdb}.pdb"
    pdb_path  = os.path.join(drc, fl_pdb)
    atoms_pdb = pr.atom.read(pdb_path)

    # Create a lookup table for this pdb...
    atom_dict = pr.atom.create_lookup_table(atoms_pdb)

    # Obtain coordinates...
    xyzs = pr.atom.extract_backbone_xyz(atom_dict, chain, nterm, cterm)

    # Calculate distance matrix...
    dmat = pr.distance.calc_dmat(xyzs, xyzs)
    dmats[i_fl, :, :] = dmat[:, :]

    ## # Export eps...
    ## dmat_bin = np.array(pr.utils.bin_image(dmat, binning = 10, nan_replace = np.nan))
    ## fl_dmat = f"{pdb}.{chain}.dmat"
    ## pal = "set palette defined ( 0 'yellow', 0 'white', 0.5 'blue', 1 'navy' )"
    ## plot_dmat(dmat_bin, fl_dmat, lbl = {}, palette = pal, smooth = True)

rmsd_dmat     = pr.distance.calc_rmsd_mats(dmats)
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

plot_dmat( rmsd_dmat_bin, 
           fl_rmsd_dmat, 
           lbl = lbls,
           lbl_font_size = 14,
           palette = pal,
           intst_max = "*",
           upper = 0.0,
           smooth = False,
           cmds_top      = cmds_guideline_top,
           cmds_bottom   = cmds_guideline_bottom,
         )

