#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' Align a list of molecules using `super` command in PyMol.  The first item 
    in the list is considered as the reference.  
'''

import pymolPy3
import pyrotein as pr
import os
from loaddata import load_gpcrdb_xlsx
## from pmview import view_dict

job_name = "xfam"

# [[[ OBTAIN THE CONSENSUS SEQUENCE ]]]
# Read the sequence alignment result...
fl_aln   = f"{job_name}.step4.psa.fil.fasta"
seq_dict = pr.fasta.read(fl_aln)

# The default ones are not good enough
nseqi = 616
cseqi = 1248

# Define atoms used for distance matrix analysis...
backbone = ["N", "CA", "C", "O"]
len_res = len(backbone)
len_seq = (cseqi - nseqi + 1) * len_res

# [[[ LOAD DATABASE ]]]
# Specify chains to process...
fl_chain = "gpcrdb.all.xlsx"
sheet    = f"{job_name}"
lines    = load_gpcrdb_xlsx(fl_chain, sheet = sheet, splitchain = True)
drc      = "pdb"

# Define atoms used for structural alignment...
backbone = ["N", "CA", "C", "O"]

# Specify the backbone atoms to select...
backbone_select = "name " + "+".join( backbone )

# Specify the rigid framework...
fl_fwk = 'fwk.dat'
fwk    = pr.utils.read_file(fl_fwk)
fwk    = [ [int(i) // 4, int(j) // 4] for i, j in fwk ]


# Pick a PDB...
pdb, chain = "6PS6", "A"
fl_pdb     = f"{pdb}.pdb"
pdb_path   = os.path.join(drc, fl_pdb)
atoms_pdb  = pr.atom.read(pdb_path)
atom_dict  = pr.atom.create_lookup_table(atoms_pdb)
chain_dict = atom_dict[chain]

# Obtain seq string for the current chain...
tar_seq = seq_dict[f"{pdb.lower()}_{chain}"]

# Obtain the mapping from seqi to resi...
seqi_to_resi_dict = pr.atom.seqi_to_resi(chain_dict, tar_seq, nseqi, cseqi)
tar_seq_fil = tar_seq[ nseqi : cseqi + 1 ]
resi_min = seqi_to_resi_dict[nseqi + pr.fasta.get_lseqi(tar_seq_fil)]
resi_max = seqi_to_resi_dict[nseqi + pr.fasta.get_rseqi(tar_seq_fil)]

# Find all resi required for alignment...
fwk_resi = []
for b, e in fwk:
    for i in range(b, e + 1):
        seqi = nseqi + i
        resi = seqi_to_resi_dict[seqi]
        if resi is None: continue
        fwk_resi.append(resi)

fwk_select = {}
fwk_select[f"{pdb}_{chain}"] = '+'.join( [ str(i) for i in fwk_resi ] )

# Customize cartoon representation...
color = { "rigid"  : "0xc1ffc1",
          "mobile" : "0xb8b8ff", }

# Start pymol
pm = pymolPy3.pymolPy3()
pm("window size, 1500, 1500")
pm("bg white")
pm("set cartoon_fancy_helices, 1")
pm("set cartoon_highlight_color, grey90")
pm("set cartoon_dumbbell_length, 1")
pm("set cartoon_dumbbell_width, 0.3")
pm("set cartoon_dumbbell_radius, 0.2")
pm("set sphere_scale, 0.3")

# Load the first structure (target)...
entry  = f"{pdb}"
pdb_path   = os.path.join(drc, f"{entry}.pdb")
pm(f"load {pdb_path}")
pm(f"remove {entry} and not chain {chain}")
pm(f"remove {entry} and not polymer.protein")
pm(f"hide cartoon, chain {chain}")
pm(f"show cartoon, chain {chain} and resi {resi_min}-{resi_max}")

# Select the rigid framework from the target...
target = f"{entry}_fwk"
pm(f"select {target}, (%{entry} and {backbone_select}) and (resi {fwk_select[f'{pdb}_{chain}']})")
pm(f"disable %{target}")

pm(f"set cartoon_color, {color['mobile']}, all")
pm(f"set cartoon_color, {color['rigid']}, {entry}_fwk")

# Create labels...
entry_dict = { f"{v[7]}_{v[10]}" : i for i, v in enumerate(lines) }

# Selectively plot entries...
pdb_list = [ "3SN6_R" ]
entry_fil_dict = { k : entry_dict[k] for k in pdb_list if k in entry_dict }

for v in pdb_list[:]:
    # Load a mobile structure...
    pdb, chain = v.split("_")
    entry      = f"{pdb}"
    pdb_path   = os.path.join(drc, f"{entry}.pdb")
    atoms_pdb  = pr.atom.read(pdb_path)
    atom_dict  = pr.atom.create_lookup_table(atoms_pdb)
    chain_dict = atom_dict[chain]

    # Obtain seq string for the current chain...
    tar_seq = seq_dict[f"{pdb.lower()}_{chain}"]

    # Obtain the mapping from seqi to resi...
    seqi_to_resi_dict = pr.atom.seqi_to_resi(chain_dict, tar_seq, nseqi, cseqi)
    tar_seq_fil = tar_seq[ nseqi : cseqi + 1 ]
    resi_min = seqi_to_resi_dict[nseqi + pr.fasta.get_lseqi(tar_seq_fil)]
    resi_max = seqi_to_resi_dict[nseqi + pr.fasta.get_rseqi(tar_seq_fil)]

    # Find all resi required for alignment...
    fwk_resi = []
    for b, e in fwk:
        for i in range(b, e + 1):
            seqi = nseqi + i
            resi = seqi_to_resi_dict[seqi]
            if resi is None: continue
            fwk_resi.append(resi)

    fwk_select[f"{pdb}_{chain}"] = '+'.join( [ str(i) for i in fwk_resi ] )

    # Work on it
    pm(f"load {pdb_path}")
    pm(f"remove {entry} and not chain {chain}")
    pm(f"remove {entry} and not polymer.protein")
    pm(f"hide cartoon, {entry} and chain {chain}")
    pm(f"show cartoon, {entry} and chain {chain} and resi {resi_min}-{resi_max}")

    # Select the rigid framework from the mobile...
    mobile = f"{entry}_fwk"
    pm(f"select {mobile}, (%{entry} and {backbone_select}) and (resi {fwk_select[f'{pdb}_{chain}']})")
    pm(f"disable %{mobile}")

    pm(f"set cartoon_color, {color['mobile']}, all")
    pm(f"set cartoon_color, {color['rigid']}, {entry}_fwk")

    # Align...
    pm(f"super {mobile}, {target}")

pm(f"orient")
input()
