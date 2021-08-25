#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pprint import pprint
from loaddata import load_xlsx
import pyrotein as pr


# [[[ OBTAIN THE CONSENSUS SEQUENCE ]]]
# Read the sequence alignment result...
# [WARNING] !!!sequence alignment is not trustworthy
fl_aln   = 'step1.all.fasta'
seq_dict = pr.fasta.read(fl_aln)

# Specify chains to process...
fl_chain = "rhodopsin.db.xlsx"
lines    = load_xlsx(fl_chain, sheet = "total", splitchain = True)
drc      = "pdb"

# Export to one single fasta that encompasses all interested chains...
fl_fasta_out = 'step2.interest.fasta'
with open(fl_fasta_out,'w') as fh:
    # Select chains according to database...
    for i_fl, line in enumerate(lines):
        # Unpack parameters
        _, pdb, chain = line[:3]

        k = f"{pdb.lower()}_{chain}"
        v = seq_dict[k]

        fh.write(f">{k}\n")
        fh.write(f"{v}\n")
