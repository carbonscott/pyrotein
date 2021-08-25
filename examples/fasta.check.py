#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from pprint import pprint
from loaddata import load_xlsx
import pyrotein as pr
import numpy as np


# [[[ OBTAIN THE CONSENSUS SEQUENCE ]]]
# Read the sequence alignment result...
# [WARNING] !!!sequence alignment is not trustworthy
fl_aln   = 'step3.msa.fasta'
seq_dict = pr.fasta.read(fl_aln)


# Obtain the consensus sequence (super seq)...
tally_dict = pr.fasta.tally_resn_in_seqs(seq_dict)
super_seq  = pr.fasta.infer_super_seq(tally_dict)


# [[[ FIND SIZE OF DISTANCE MATRIX ]]]
# Get the sequence index (alignment) on the n-term side...
nseqi = pr.fasta.get_lseqi(super_seq)
cseqi = pr.fasta.get_rseqi(super_seq)
## cseqi = 290

# Define atoms used for distance matrix analysis...
backbone = ["N", "CA", "C", "O"]
len_res = len(backbone)
len_seq = (cseqi - nseqi + 1) * len_res

# [[[ RMSD ANALYSIS ]]]
drc        = "pdb"
pdb, chain = "5dgy", "C"
fl_pdb     = f"{pdb}.pdb"
pdb_path   = os.path.join(drc, fl_pdb)
atoms_pdb  = pr.atom.read(pdb_path)
atom_dict  = pr.atom.create_lookup_table(atoms_pdb)
chain_dict = atom_dict[chain]

# Obtain seq string for the current chain...
tar_seq = seq_dict[f"{pdb.lower()}_{chain}"]

# Obtain xyzs
entry = f"{pdb}_{chain}"
print(entry)
## xyzs = pr.atom.extract_xyz_by_seq(backbone, chain_dict, tar_seq, nseqi, cseqi)


if True:
    # Extract resn to resi mapping (like the sequence viewer on PyMol)...
    # Non amino acid residue (like ligand) are bypassed
    resn_to_resi_dict = pr.atom.resn_to_resi(chain_dict)

    # Select the bound sequence by nseqi and cseqi...
    tar_seq_bound           = tar_seq[nseqi : cseqi + 1]
    tar_seq_bound_continous = tar_seq[nseqi : cseqi + 1].replace('-', '')

    # Obtain the original sequence from PDB...
    # May have overhead in the beginning or the end
    seq_orig = ''.join([ v for v in resn_to_resi_dict.values() ])

    # Obtain the starting index by string match...
    lb_term = seq_orig.find(tar_seq_bound_continous)

    # Obtain the ending index by the length of the coutinous (no '-') sequence...
    ub_term = lb_term + len(tar_seq_bound_continous)

    # Obtain list of resi bound by nseqi and cseqi...
    resi_list       = [ v for v in resn_to_resi_dict.keys() ]
    resi_bound_list = resi_list[lb_term : ub_term]

    # Initialize mapping...
    seqi_to_resi_dict = { k : None for k in range(nseqi, cseqi + 1) }

    # Counter to go through the bound sequence by nseqi and cseqi...
    res_counter = 0

    # Loop through
    for i, seqi in enumerate(range(nseqi, cseqi + 1)):
        # Skip the '-' residue...
        if tar_seq_bound[i] == '-': continue

        # Access the resi...
        resi = resi_bound_list[res_counter]

        # Record the mapping...
        seqi_to_resi_dict[seqi] = resi

        # Increment the residue counter...
        res_counter += 1


