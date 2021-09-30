#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' Very ad-hoc code. Probably won't make it any better.
'''

import pymolPy3
import pyrotein as pr
import os
from loaddata import load_gpcrdb_xlsx
import numpy as np
import givens as gv
import colorsimple as cs
from display import calc_rigid_framework

# [[[ User interface ]]]
job_name = "xfam-loop"
pdb1 = "6X18_R".split("_")
pdb2 = "6X19_R".split("_")
clr1 = "0x00A600"
clr2 = "red"
pdb1[0] = pdb1[0].lower()
pdb2[0] = pdb2[0].lower()


# [[[ Hide from users ]]]

# Read the sequence alignment result...
fl_aln   = f"xfam.step4.psa.fil.fasta"
seq_dict = pr.fasta.read(fl_aln)

# Obtain the consensus sequence (super seq)...
tally_dict = pr.fasta.tally_resn_in_seqs(seq_dict)
super_seq  = pr.fasta.infer_super_seq(tally_dict)

# Specify chains to process...
fl_chain = "gpcrdb.all.xlsx"
sheet    = f"xfam"
lines    = load_gpcrdb_xlsx(fl_chain, sheet = sheet, splitchain = True)
drc      = "pdb"
pdb1_check = pdb1[0].upper()
pdb2_check = pdb2[0].upper()
lines = list(filter( lambda x: x[7] in [pdb1_check, pdb2_check], lines ))


# Correct order in lines...
if pdb1 != lines[0][7]: lines = lines[::-1]

# The default ones are not good enough
nseqi = 637
cseqi = 1371
num_seq = cseqi - nseqi + 1

# Define atoms used for distance matrix analysis...
backbone = ["N", "CA", "C", "O"]
len_res = len(backbone)
len_seq = (cseqi - nseqi + 1) * len_res

# Accumulate distance matices as lower triangluar matrix...
len_lower_tri = (len_seq * len_seq - len_seq) // 2
dmats = np.zeros((len(lines), len_lower_tri))

# Build a mask to filter out loop segments...
mask_loop = np.zeros((len_seq, len("xyz")), dtype = np.int8)

# Use labels to create mask...
# Notice the +1 for right-inclusion is removed
labels = {'H8': [722, 735],
          'TM1': [0, 34],
          'TM2': [63, 96],
          'TM3': [158, 194],
          'TM4': [227, 254],
          'TM5': [316, 357],
          'TM6': [591, 631],
          'TM7': [680, 708]}
for k, v in labels.items(): labels[k] = [ (i) * len_res for i in v ]

# Create mask...
# +1 offset is to deal with numpy indexing which is not right-inclusive 
for _, (b, e) in labels.items(): mask_loop[b:e+1, :] = 1

drc_dmat = "xfam-loop.dmats"
pal = "set palette defined \
       (-10 '#800000', -5 'red', -1 'white', 0 'seagreen', \
          1 'white'  , 5 'blue', 10 'navy')"

for i_fl, line in enumerate(lines[:]):
    # Unpack parameters
    uniprot = line[1].lower()
    spec    = line[5].lower()
    pdb     = line[7].lower()
    chain   = line[10].upper()

    # Read coordinates from a PDB file...
    fl_pdb    = f"{pdb}.pdb"
    pdb_path  = os.path.join(drc, fl_pdb)
    atoms_pdb = pr.atom.read(pdb_path)

    # Create a lookup table for this pdb...
    atom_dict = pr.atom.create_lookup_table(atoms_pdb)

    # Obtain the chain to process...
    chain_dict = atom_dict[chain]

    # Obtain seq string for the current chain...
    tar_seq = seq_dict[f"{pdb.lower()}_{chain}"]

    # Obtain xyzs...
    entry = f"{pdb}_{chain}"
    xyzs = pr.atom.extract_xyz_by_seq(backbone, chain_dict, tar_seq, nseqi, cseqi)
    xyzs[mask_loop == 0] = np.nan

    # Calculate distance matrix...
    dmat = pr.distance.calc_dmat(xyzs, xyzs)

    # Calculate distance matrix...
    dmat = pr.distance.calc_dmat(xyzs, xyzs)
    dmats[i_fl, :] = pr.utils.mat2tril(dmat, offset = -1)

# Replace np.nan with mean across samples...
pr.utils.fill_nan_with_mean(dmats.T, axis = 1)
pr.utils.fill_nan_with_zero(dmats)
dmat_diff = pr.distance.calc_rmsd_array(dmats)

# Rescale (normalize) by the minmax...
minv, maxv = np.min(dmat_diff), np.max(dmat_diff)
dmat_diff = (dmat_diff - minv) / (maxv - minv)

# Convert `u` at position `rank` into a lower triangular matrix...
dmat_diff_2d = pr.utils.array2tril(dmat_diff, len_seq, offset = -1)

# Restore dmat to a full matrix for visualization...
dmat_diff_2d = dmat_diff_2d + dmat_diff_2d.T

# Read the initial framework...
# +1 to right-include the index has been considered in the fl_fwk
fl_fwk = f"fwk.{job_name}.dat"
seqi_fwk_list = pr.utils.read_file(fl_fwk, str_to_num = True, num_type = int)
fwk_list = calc_rigid_framework(dmat_diff_2d, seqi_fwk_list, num_seq, len_res, min_size = 5, min_mean = 0.001, epsilon = 0.0001)

fwk = []
for i in range(len(fwk_list)):
    b1, e1 = [ fwk_list[i][0] // len_res, fwk_list[i][-1] // len_res ]
    fwk.append((b1, e1))


# [[[ LOAD DATABASE ]]]
# Specify the backbone atoms to select...
backbone_select = "name " + "+".join( backbone )


# ======== MOLECLUE REFERENCE =========
# Pick a PDB...
## postfix = "_GPCRDB"
postfix = ""
pdb, chain = pdb1
entry = f"{pdb}{postfix}"
## drc = "pdb.numbered"
drc = "pdb"
fl_pdb     = f"{entry}.pdb"
pdb_path   = os.path.join(drc, fl_pdb)
atoms_pdb  = pr.atom.read(pdb_path)
atom_dict  = pr.atom.create_lookup_table(atoms_pdb)
chain_dict = atom_dict[chain]

# Obtain seq string for the current chain...
tar_seq = seq_dict[f"{pdb}_{chain}"]

# Obtain the mapping from seqi to resi...
seqi_to_resi_dict = pr.atom.seqi_to_resi(chain_dict, tar_seq, nseqi, cseqi)
tar_seq_fil = tar_seq[ nseqi : cseqi + 1 ]
resi_min = seqi_to_resi_dict[nseqi + pr.fasta.get_lseqi(tar_seq_fil)]
resi_max = seqi_to_resi_dict[nseqi + pr.fasta.get_rseqi(tar_seq_fil)]

# Find all resi required for alignment...
fwk_resi = []
for b, e in fwk:
    for i in range(b, e + 1):
        seqi_fwk = nseqi + i
        resi_fwk = seqi_to_resi_dict[seqi_fwk]
        if resi_fwk is None: continue
        fwk_resi.append(resi_fwk)
fwk_select = {}
fwk_select[f"{pdb}_{chain}"] = '+'.join( [ str(i) for i in fwk_resi ] )


# Start pymol
pm = pymolPy3.pymolPy3()
pm("bg white")
pm("set cartoon_fancy_helices, 1")
pm("set cartoon_highlight_color, grey90")
pm("set cartoon_dumbbell_length, 1")
pm("set cartoon_dumbbell_width, 0.3")
pm("set cartoon_dumbbell_radius, 0.2")
pm("set sphere_scale, 0.3")

# Load the first structure (target)...
pdb_path   = os.path.join(drc, f"{entry}.pdb")
pm(f"load {pdb_path}")
pm(f"set cartoon_color, {clr1}, %{entry}")
pm(f"remove {entry} and not chain {chain}")
pm(f"remove {entry} and not polymer.protein")
pm(f"hide cartoon, chain {chain}")
pm(f"show cartoon, chain {chain} and resi {resi_min}-{resi_max}")

# Select the rigid framework from the target...
target = f"{entry}_fwk"
pm(f"select {target}, (%{entry} and {backbone_select}) and (resi {fwk_select[f'{pdb}_{chain}']})")
pm(f"set cartoon_color, tv_yellow, %{target}")
pm(f"disable %{target}")


# ======== MOLECLUE MOBILE 1 =========
pdb, chain = pdb2
entry      = f"{pdb}{postfix}"
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
        seqi_fwk = nseqi + i
        resi_fwk = seqi_to_resi_dict[seqi_fwk]
        if resi_fwk is None: continue
        fwk_resi.append(resi_fwk)

fwk_select[f"{pdb}_{chain}"] = '+'.join( [ str(i) for i in fwk_resi ] )

# Work on it
pm(f"load {pdb_path}")
pm(f"set cartoon_color, {clr2}, %{entry}")
pm(f"remove {entry} and not chain {chain}")
pm(f"remove {entry} and not polymer.protein")
pm(f"hide cartoon, {entry} and chain {chain}")
pm(f"show cartoon, {entry} and chain {chain} and resi {resi_min}-{resi_max}")

# Select the rigid framework from the mobile...
mobile = f"{entry}_fwk"
pm(f"select {mobile}, (%{entry} and {backbone_select}) and (resi {fwk_select[f'{pdb}_{chain}']})")
pm(f"set cartoon_color, tv_yellow, %{mobile}")
pm(f"disable %{mobile}")

# Align...
pm(f"super {mobile}, {target}")


pm("orient")
pm("window size, 1500, 1500")

input()
