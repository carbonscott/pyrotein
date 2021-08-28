#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pyrotein as pr
from loaddata import load_gpcrdb_xlsx
import pymolPy3

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
    drc = "pdb.refine"
    fl_dat = "pdb.refine.dat"
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
