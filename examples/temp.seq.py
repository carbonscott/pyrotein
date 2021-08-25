#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## import sys
## sys.path.insert(0, '/home/scott/Dropbox/codes/pyrotein')

import pyrotein as pr
import os

fl_aln = 'seq.align.fasta'
seq_dict = pr.fasta.read(fl_aln)

tally_dict = pr.fasta.tally_resn_in_seqs(seq_dict)
super_seq = pr.fasta.infer_super_seq(tally_dict)
seq_to_resi_dict = pr.fasta.seq_to_resi(super_seq, 1)

ref = super_seq
pdb = '1f88'
chain = 'A'
entry = f"{pdb}_{chain}"
tar = seq_dict[entry]

seq_diff = pr.fasta.diff_seq(tar, ref)

nterm, cterm = 1, 322

ref_simp = pr.fasta.strip_null(ref)
seq_to_resi_dict = pr.fasta.seq_to_resi(ref_simp, 1)
nseqi    = pr.fasta.get_lseqi(tar)
cseqi    = pr.fasta.get_rseqi(tar)
tar_simp = tar[nseqi : nseqi + len(ref_simp)]
seq_simp_dict = pr.fasta.seq_to_resi(ref_simp, 1)
seq_simp_diff = pr.fasta.diff_seq(tar_simp, ref_simp)
seq_non_null_list = pr.fasta.seqi_non_null(seq_simp_diff)


# Read coordinates from a PDB file...
fl_pdb    = f"{pdb}.pdb"
drc       = 'pdb'
pdb_path  = os.path.join(drc, fl_pdb)
atoms_pdb = pr.atom.read(pdb_path)

# Create a lookup table for this pdb...
atom_dict = pr.atom.create_lookup_table(atoms_pdb)

