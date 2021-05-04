#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## import sys
## sys.path.insert(0, '/home/scott/Dropbox/codes/pyrotein')

import numpy as np
import pyrotein as pr
import os
from loaddata import load_xlsx
from display import showHistogram

if True:
    # [[[ OBTAIN THE CONSENSUS SEQUENCE ]]]
    # Read the sequence alignment result...
    # [WARNING] !!!sequence alignment is not trustworthy, need to check manually
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

    # Create a dmask...
    dmask = pr.utils.sparse_mask(super_seg)

    # Load constant -- atomlabel...
    label_dict = pr.atom.constant_atomlabel()
    aa_dict    = pr.atom.constant_aminoacid_code()

    # Calculate the total length of distance matrix...
    len_dmat = np.sum( [ len(label_dict[aa_dict[i]]) for i in super_seg ] )

    # [[[ ANALYZE PDB ENTRIES ]]]
    # Specify chains to process...
    fl_chain = "chains.comp.xlsx"
    lines    = load_xlsx(fl_chain, sheet = "Sheet1")
    drc      = "pdb"
    len_lower_tri = (len_dmat * len_dmat - len_dmat) // 2
    dmats = np.zeros((len(lines), len_lower_tri))

    # Create sparse mask...
    dmask = pr.utils.sparse_mask(super_seg, offset = 0) * \
            pr.utils.sparse_mask(super_seg, offset = 1, val = 0.0)

    # Process each entry...
    for i_fl, line in enumerate(lines[:]):
        # Unpack parameters
        _, pdb, chain, species = line[:4]

        print(f"Processing {pdb}_{chain}")

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

        # Apply sparse matrix mask...
        dmat *= dmask

        # Convert dmat into one-dimensional array and keep it in dmats...
        dmats[i_fl, :] = pr.utils.mat2tril(dmat, offset = -1)

    # [DO IT LATER] Replace np.nan with mean across samples...
    pr.utils.fill_nan_with_zero(dmats)

    # Set cutoff range...
    dmin, dmax = 1.0, 5.0

    # Find the indices of values that are within the range of (dmin, dmax)...
    in_range_bool      = np.logical_and( dmin < dmats, dmats < dmax )
    in_range_bool_flat = np.any(in_range_bool, axis = 0)
    in_range_index     = np.argwhere( in_range_bool_flat ).reshape(-1)

    # Form a submatrix by selecting values within the range...
    dmats_sub = dmats[:, in_range_index]

## dmats_sub = np.load("dmats_sub.npy")

cmds = [ "set xtics 0.1 rotate by 90 right font ',10'" ]
showHistogram(dmats_sub, bin_cap = 700, rng = [dmin, dmax], title = "contact", cmds = cmds)
