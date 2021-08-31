#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from pprint import pprint
import pyrotein as pr
import numpy as np


# [[[ OBTAIN THE CONSENSUS SEQUENCE ]]]
# Read the sequence alignment result...
## fl_aln   = 'step3.psa.fasta'
fl_aln   = 'xfam.step4.psa.fil.fasta'
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
pdb, chain = "3cap", "A"
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


def get_seg_nterm(seq, seg): return seq.find(seg)
def get_seg_cterm(seq, seg): return seq.find(seg) + len(seg) - 1

TM_dict = {
    "TM1" : { "nterm" : "EPWQ", "cterm" : "TVQH", },
    "TM2" : { "nterm" : "TPLN", "cterm" : "TSLHG", },
    "TM3" : { "nterm" : "GPTGCN", "cterm" : "VVVCK", },
    "TM4" : { "nterm" : "GENHA", "cterm" : "PPLVG", },
    "TM5" : { "nterm" : "NNESF", "cterm" : "AAAQQ", },
    "TM6" : { "nterm" : "SATTQK", "cterm" : "FYIFTH", },
    "TM7" : { "nterm" : "GPIFMT", "cterm" : "YIMM", },
    "H8"  : { "nterm" : "NKQF", "cterm" : "MVTTLC", },
}

nseqi = get_seg_nterm(tar_seq, TM_dict["TM1"]["nterm"])
cseqi = get_seg_cterm(tar_seq, TM_dict["H8" ]["cterm"])

labels = { }
for k, v in TM_dict.items():
    labels[k] = [ 
        get_seg_nterm(tar_seq, v["nterm"]) - nseqi,
        get_seg_cterm(tar_seq, v["cterm"]) + 1 - nseqi,
    ]


pprint(labels)
