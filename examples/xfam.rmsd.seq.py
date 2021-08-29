#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' Replace 'X' in fl_aln with '-'
'''

## import sys
## sys.path.insert(0, "/home/scott/Dropbox/codes/pyrotein")
## sys.path.insert(0, "/Users/scott/Dropbox/codes/pyrotein")

import os
import numpy as np
import pyrotein as pr
from loaddata import load_gpcrdb_xlsx
from display import plot_dmat

# Set job name...
job_name = "xfam"

# [[[ OBTAIN THE CONSENSUS SEQUENCE ]]]
# Read the sequence alignment result...
fl_aln   = f"{job_name}.step4.psa.fil.fasta"
seq_dict = pr.fasta.read(fl_aln)

# Obtain the consensus sequence (super seq)...
tally_dict = pr.fasta.tally_resn_in_seqs(seq_dict)
super_seq  = pr.fasta.infer_super_seq(tally_dict)

# [[[ LOAD DATABASE ]]]
# Specify chains to process...
fl_chain = "gpcrdb.all.xlsx"
sheet    = f"{job_name}"
lines    = load_gpcrdb_xlsx(fl_chain, sheet = sheet, splitchain = True)
drc      = "pdb"

# [[[ FIND SIZE OF DISTANCE MATRIX ]]]
# Get the sequence index (alignment) on the n-term side...
## nseqi = pr.fasta.get_lseqi(super_seq)
## cseqi = pr.fasta.get_rseqi(super_seq)

# The default ones are not good enough
nseqi = 616
cseqi = 1248

# Define atoms used for distance matrix analysis...
backbone = ["N", "CA", "C", "O"]
len_res = len(backbone)
len_seq = (cseqi - nseqi + 1) * len_res

# Preallocate memory for dmats...
dmats = np.zeros((len(lines), len_seq, len_seq))

# [[[ PROCESS RMSD ]]]
for i_fl, line in enumerate(lines):
    # Unpack parameters
    uniprot = line[1].lower()
    spec    = line[5].lower()
    pdb     = line[7].lower()
    chain   = line[10].upper()

    # Read coordinates from a PDB file...
    fl_pdb    = f"{pdb}.pdb"
    pdb_path  = os.path.join(drc, fl_pdb)
    atoms_pdb = pr.atom.read(pdb_path)

    # Create a lookup table for this pdb...
    atom_dict = pr.atom.create_lookup_table(atoms_pdb)

    # Obtain the chain to process...
    chain_dict = atom_dict[chain]

    # Obtain seq string for the current chain...
    tar_seq = seq_dict[f"{pdb.lower()}_{chain}"]

    # Obtain xyzs
    entry = f"{pdb}_{chain}"
    print(entry)
    xyzs = pr.atom.extract_xyz_by_seq(backbone, chain_dict, tar_seq, nseqi, cseqi)

    # Calculate distance matrix...
    dmat = pr.distance.calc_dmat(xyzs, xyzs)
    dmats[i_fl, :, :] = dmat[:, :]

rmsd_dmat = pr.distance.calc_rmsd_mats(dmats)
np.save("rmsd_dmat.seq.npy" , rmsd_dmat)
np.save("rmsd_len.seq.npy"  , len(backbone))
np.save("nseqi.seq.npy", nseqi)
np.save("cseqi.seq.npy", cseqi)
