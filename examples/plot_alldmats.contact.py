#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, "/home/scott/Dropbox/codes/pyrotein")
## sys.path.insert(0, "/Users/scott/Dropbox/codes/pyrotein")

import os
import numpy as np
import pyrotein as pr
from loaddata import load_xlsx, label_TMs
from display import plot_dmat
import multiprocessing as mp

# [[[ OBTAIN THE CONSENSUS SEQUENCE ]]]
# Read the sequence alignment result...
# [WARNING] !!!sequence alignment is not trustworthy
fl_aln   = 'seq.align.fasta'
seq_dict = pr.fasta.read(fl_aln)

# Obtain the consensus sequence (super seq)...
tally_dict = pr.fasta.tally_resn_in_seqs(seq_dict)
super_seq  = pr.fasta.infer_super_seq(tally_dict)


# [[[ FIND SIZE OF DISTANCE MATRIX ]]]
# Get the sequence index (alignment) on the n-term side...
nseqi = pr.fasta.get_lseqi(super_seq)

# User defined range...
nterm, cterm = 1, 322
len_seg = cterm - nterm + 1
super_seg = super_seq[nseqi : nseqi + len_seg]

# Create a dmat mask...
dmask = pr.utils.sparse_mask(super_seg)

# [[[ ANALYZE PDB ENTRIES ]]]
# Specify chains to process...
fl_chain = "chains.comp.xlsx"
lines    = load_xlsx(fl_chain, sheet = "Sheet1")
drc      = "pdb"
drc_dmat = "dmats.full"

# Create labels for u...
labels_TM = label_TMs()
for k, v in labels_TM.items(): 
    labels_TM[k] = [ pr.fasta.resi_to_seqi(i, super_seg, nterm) for i in v ]

## pal = "set palette defined ( 0 '#F6FF9E', 0 'white', 0.5 'blue', 1 'navy' )"
## pal = '''
## set palette negative defined ( \
##     0 '#D53E4F',\
##     1 '#F46D43',\
##     2 '#FDAE61',\
##     3 '#FEE08B',\
##     4 '#E6F598',\
##     5 '#ABDDA4',\
##     6 '#66C2A5',\
##     7 '#3288BD' )
## '''
pal = '''
set palette defined ( 0 "#000090",\
                      1 "#000fff",\
                      2 "#0090ff",\
                      3 "#0fffee",\
                      4 "#90ff70",\
                      5 "#ffee00",\
                      6 "#ff7000",\
                      7 "#ee0000",\
                      8 "#7f0000")
'''

for i_fl, line in enumerate(lines[-1:]):
    # Unpack parameters
    _, pdb, chain, _ = line[:4]

    # Read coordinates from a PDB file...
    fl_pdb    = f"{pdb}.pdb"
    pdb_path  = os.path.join(drc, fl_pdb)
    atoms_pdb = pr.atom.read(pdb_path)

    # Create a lookup table for this pdb...
    atom_dict = pr.atom.create_lookup_table(atoms_pdb)

    # Obtain the target protein by range...
    tar_seq = seq_dict[f"{pdb}_{chain}"]
    tar_seg = tar_seq[nseqi : nseqi + len_seg]

    # Standardize sidechain atoms...
    pr.atom.standardize_sidechain(atom_dict)

    # Obtain coordinates...
    xyzs = pr.atom.extract_xyz_by_seq(tar_seg, super_seg, atom_dict, chain, nterm, cterm)

    # Calculate distance matrix...
    dmat = pr.distance.calc_dmat(xyzs, xyzs)
    pr.utils.fill_nan_with_zero(dmat)

    # Set cutoff range...
    dmin, dmax = 2.65, 3.0

    # Find the indices of values that are within the range of (dmin, dmax)...
    out_range_bool      = np.logical_or( dmin > dmat, dmat > dmax )

    # Form a submatrix by selecting values within the range...
    dmat[out_range_bool] = np.nan

    ## # Apply dmask...
    ## dmat *= dmask

    ## # Manual intervention (trivial)...
    ## fl_exclude = os.path.join(drc_dmat, f"{pdb}.{chain}.dmat.exclude.dat")
    ## exclude_list = pr.utils.read_file(fl_exclude, numerical = True)
    ## exclude_list = [ [int(x), int(y)] for x, y in exclude_list ]
    ## for x, y in exclude_list: dmat[x,y] = np.nan

    ## # Manual intervention (H bond)...
    ## fl_exclude = os.path.join(drc_dmat, f"{pdb}.{chain}.dmat.exclude.nonHbond.dat")
    ## exclude_list = pr.utils.read_file(fl_exclude, numerical = True)
    ## exclude_list = [ [int(x), int(y)] for x, y in exclude_list ]
    ## for x, y in exclude_list: dmat[x,y] = np.nan

    # Put more resi labels...
    diaglbl = {}
    for resi in range(nterm, cterm, 1):
        seqi = pr.fasta.resi_to_seqi(resi, super_seg, nterm)
        diaglbl[resi] = [ seqi, seqi ]

    fl_dmat = os.path.join(drc_dmat, f"{pdb}.{chain}.dmat")
    plot_dmat(dmat, 
              fl_dmat, 
              lbl = labels_TM,
              lbl_fontsize = 18,
              diaglbl = diaglbl,
              diaglblfontsize = 3,
              intst_min = dmin,
              intst_max = dmax,
              palette = pal, 
              temp    = False,
              mode    = 'sparse',
              showsparselabel = False,
              NaN = 0)
