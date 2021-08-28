#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pyrotein as pr

drc_i = "psa.i"
drc_o = "psa.o"
fl_dat = "v1.s4.req.refine.dat"
lines = pr.utils.read_file(fl_dat)

def fli_for_needle(pdb, chain):
    return [
        f"tmp.{pdb}_{chain}.fasta",
    ]

def flo_from_needle(pdb, chain):
    return [
        f"out.{pdb}_{chain}.aln.txt",
        f"out.{pdb}_{chain}.bsequence.txt",
        f"out.{pdb}_{chain}.asequence.txt",
        f"out.{pdb}_{chain}.submission.params",
    ]

for line in lines:
    pdb, chain = line
    for fl in fli_for_needle(pdb, chain): 
        path_fl = os.path.join(drc_i, fl)
        print(f"Removing {path_fl}")
        if os.path.exists(path_fl): os.remove(path_fl)
    for fl in flo_from_needle(pdb, chain): 
        path_fl = os.path.join(drc_o, fl)
        print(f"Removing {path_fl}")
        if os.path.exists(path_fl): os.remove(path_fl)
