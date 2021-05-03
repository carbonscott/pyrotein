#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## import sys
## sys.path.insert(0, "/home/scott/Dropbox/codes/pyrotein")
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

# Build a list of (resi, resn, atom)...
label_list = pr.utils.label_dmat(super_seg, nterm, cterm)


# Obtain the dist_list...
pdb = '6qno'
chain = 'R'
drc = "dmats.full"
dmin, dmax = 2.65, 3.0
fl_u = os.path.join(drc, f"{pdb}.{chain}.dmat.dat")
dist_list = pr.utils.read_file(f"{fl_u}", numerical = True)

def filter_dist(d_list, label_list, dmin, dmax):
    d_list_filter = [ (int(x), int(y)) for x,y,v in d_list if dmin < v < dmax ]
    d_set = set([ j for i in d_list_filter for j in i ])

    d_list_atm = [ (label_list[int(x)], label_list[int(y)], v) for x,y,v in d_list if dmin < v < dmax ]
    return sorted(list(d_set)), sorted(d_list_atm, key = lambda x: x[2])

dist_set, dist_pair = filter_dist(dist_list, label_list, dmin, dmax)

for x,y,v in dist_pair: print(f"{x:14s}  {y:14s}  {v:.4f}")

atm_pair = []
for x,y,v in dist_pair:
    id_x, _, atm_x = x.split('.')
    id_y, _, atm_y = y.split('.')
    atm_pair.append(f"{atm_x}-{atm_y}:{int(id_y) - int(id_x):1d}")

tally_atm_pair = pr.utils.tally_int(atm_pair)
tally_atm_pair_sorted = {k: v for k, v in sorted(tally_atm_pair.items(), key=lambda item: item[1], reverse = True)}

def print_tally(tally_atm_pair_sorted, top = 5):
    print("---- Tally ----")
    for i, (k, v) in enumerate(tally_atm_pair_sorted.items()):
        if i >= top: break
        print(f"{k:8s} => {v:6d}")

print_tally(tally_atm_pair_sorted, top = 5)


def print_pair(dist_list, atm1, atm2):
    for x,y,v in dist_pair:
        id_x, _, atm_x = x.split('.')
        id_y, _, atm_y = y.split('.')
        if atm_x != atm1 or atm_y != atm2: continue
        print(f"{x:14s} {y:14s}  {v:.4f}")
