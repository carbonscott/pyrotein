#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pyrotein as pr
from loaddata import load_gpcrdb_xlsx, load_uniprot_species, name_to_code

# Read uniprot specification...
fl_spec   = "organism.dat"
spec_dict = load_uniprot_species(fl_spec)
name_to_code_dict = name_to_code(spec_dict, "c")    # "c" stands for common name

# Specify chains to process...
fl_xlsx    = "gpcrdb.all.xlsx"
sheet      = "total"
lines      = load_gpcrdb_xlsx(fl_xlsx, sheet = sheet)
drc        = "pdb"
drc_input  = "psa.i"
drc_export = "psa.o"

for line in lines[:]:
    # Unpack parameters
    uniprot = line[1].lower()
    spec    = line[5].lower()
    pdb     = line[7].lower()
    chain   = line[10].upper()

    # Form entry name...
    entry = f"{pdb}_{chain}"

    # Derive the id like uniprot_speccode...
    if not spec in name_to_code_dict: 
        print(f"!!! UniProt entry error -- No {spec}")
        continue
    speccode = name_to_code_dict[f"{spec}"]
    entry_ref = f"{uniprot}_{speccode}"
    print(f"Processing {entry} from {entry_ref}...")

    fl_export = os.path.join(drc_export, f"out.{entry}.aln.txt")
    if os.path.exists(fl_export): seq_dict = pr.fasta.read(fl_export)
    else: 
        print(f"!!!{fl_export} not exists")
        continue

    if not entry_ref in seq_dict:
        print(f"{entry_ref} N/A!!!  Skip this species...")
        continue

    seq_disp = seq_dict[entry_ref]
    if seq_disp.find("-") > 0: 
        print(f"??? {entry} having unconventional sequence")
        print(seq_disp)
