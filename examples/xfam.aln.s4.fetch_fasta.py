#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' cp drc/*.pdb   pdb/
    cp drc/*.fasta pdb/

    cat pdb/????.fasta > step1.all.fasta
    ./fasta.step1.py

    ./aln.s5.remove_defect_fasta.py
    ./aln.s1.psa.needle.py

'''

import pyrotein as pr
from loaddata import load_gpcrdb_xlsx
import pymolPy3

# Give a name to the analysis...
job_name = "xfam"

# Launch pymol without GUI...
pm = pymolPy3.pymolPy3(0)

if False:
    # Specify chains to process...
    fl_chain = "gpcrdb.all.xlsx"
    lines    = load_gpcrdb_xlsx(fl_chain, sheet = "total", splitchain = True)
    drc      = "pdb"

    # Flatten the list
    entries = [ i[7] for i in lines ]

else:
    drc = f"{job_name}.pdb.refine"
    fl_dat = f"{job_name}.req_refine.dat"
    lines = pr.utils.read_file(fl_dat)
    entries = [ i[0] for i in lines ]


# Fetching...
pm(f"cd {drc}")
for i in entries:
    pdb = i.lower()
    pm(f"load {pdb}.pdb")
    pm(f"save {pdb}.fasta")
    pm(f"delete {pdb}")
pm("quit")
