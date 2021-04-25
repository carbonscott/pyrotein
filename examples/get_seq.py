#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from loaddata import load_xlsx
import pyrotein as pr
import os


# Specify chains to process...
fl_chain = "chains.comp.xlsx"
lines    = load_xlsx(fl_chain)
drc      = "pdb"

seq_acc = []
for i_fl, line in enumerate(lines):
    # Unpack parameters
    _, pdb, chain, species = line[:4]

    # Extract sequence...
    fl_fasta = os.path.join(drc, f"{pdb}.fasta")
    seq = pr.fasta.read(fl_fasta)

    # Form a fasta string in response to chain...
    # The format only works for fasta converted from pdb by pymol
    k = f"{pdb}_{chain}"
    seq_acc.append(f">{k}")
    seq_acc.append( seq[k] )

fl_export = "seq.fasta"
with open(fl_export,'w') as fh: fh.write("\n".join(seq_acc))
