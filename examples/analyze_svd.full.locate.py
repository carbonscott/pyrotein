#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## import sys
## sys.path.insert(0, "/home/scott/Dropbox/codes/pyrotein")
## sys.path.insert(0, "/Users/scott/Dropbox/codes/pyrotein")

import os
import numpy as np
import pyrotein as pr
import GnuplotPy3

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

# [[[ ANALYZE SVD COMPONENTS ]]]
# Specify the component...
comp     = "u02"
domain   = "neg"
drc = "svd.full.u"
fl_comp = os.path.join(drc, f"{comp}.dat")

# Read the interesting dots in u plot...
dist_list = pr.utils.read_file(f"{fl_comp}", numerical = True)

# Tally residues of interaction...
# Calling it pm is because of the visualization of it in pymol
# d_list to support sorting by distance...
pm_list = []
d_list  = []
for (x, y, d) in dist_list:
    x, y = int(x), int(y)
    resi1, resn1, atom1 = label_list[int(x)].split('.')
    resi2, resn2, atom2 = label_list[int(y)].split('.')
    pm_list.extend((int(resi1), int(resi2)))
    d_list.append((int(resi1), int(resi2), d))

# Sort d_list by distance...
d_list.sort( key = lambda x: x[2] )
fl_d = f'locate.dist.{comp}.{domain}.dat'
with open(fl_d,'w') as fh:
    for x, y, d in d_list:
        fh.write(f"{x:4d}  {y:4d}  {d:.4f}")
        fh.write("\n")

# Tally...
tmin, tmax = 200, 10000
pm_filter_list = filter( lambda x: tmin < x < tmax, pm_list )
pm_tally_dict = pr.utils.tally_int(pm_filter_list)


# Visualize the tally...
gp = GnuplotPy3.GnuplotPy3()
gp("set terminal postscript eps  size 3.5, 2.62 \\")
gp("                             enhanced color \\")
gp("                             font 'Helvetica,12' \\")
gp("                             linewidth 1")
gp(f"set output 'locate.tally.{comp}.{domain}.eps'")
gp(f"set log y")

gp(f"plot \\")
gp(f"'-' using 1:2 with points pointtype 6 linecolor rgb 'black' notitle, \\")
gp(f"'-' using 1:2:3 with labels offset 2.0,0 font ',10' notitle, \\")
gp(f"")

for k, v in pm_tally_dict.items():
    gp(f"{k} {v}")
gp("e")
for k, v in pm_tally_dict.items():
    gp(f"{k} {v} {k}")
gp("e")
gp("exit")
