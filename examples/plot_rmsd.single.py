#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pyrotein as pr
from loaddata import load_xlsx, label_TMs
from display import plot_dmat, plot_simple_dmat


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
for i_fl, line in enumerate(lines[-6:-1]):
    # Unpack parameters
    _, pdb, chain, species = line[:4]

    # Read coordinates from a PDB file...
    fl_pdb    = f"{pdb}.pdb"
    pdb_path  = os.path.join(drc, fl_pdb)
    atoms_pdb = pr.atom.read(pdb_path)

    # Create a lookup table for this pdb...
    atom_dict = pr.atom.create_lookup_table(atoms_pdb)

    # Obtain coordinates...
    xyzs = pr.atom.extract_xyz(peptide, atom_dict, chain, nterm, cterm)

    # Calculate distance matrix...
    dmat = pr.distance.calc_dmat(xyzs, xyzs)
    dmats[i_fl, :, :] = dmat[:, :]

    ## # Export eps...
    ## dmat_bin = np.array(pr.utils.bin_image(dmat, binning = 10, nan_replace = np.nan))
    fl_dmat = f"{pdb}.{chain}.dmat"
    pal = "set palette defined ( 0 'yellow', 0 'white', 0.5 'blue', 1 'navy' )"
    ## plot_dmat(dmat_bin, fl_dmat, lbl = {}, palette = pal, smooth = True)
    plot_simple_dmat( dmat, 
           fl_dmat, 
           width = 12,
           height = 12,
           fontsize = 29,
           lbl = {},
           lbl_fontsize = 29,
           linewidth = 2.0,
           palette = pal,
           intst_max = "*",
           upper = np.nanmax(dmat) * 0.0125,
           smooth = False,
         )


## rmsd_dmat = pr.distance.calc_rmsd_mats(dmats)
## np.save("rmsd_dmat.npy", rmsd_dmat)
