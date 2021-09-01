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


def reverse_sign(u, vh, rank, index_from_zero = True):
    # Comply with the convention (1-based index)
    rank_in_data = rank if index_from_zero else rank - 1

    # Reverse sign...
    u[:, rank_in_data]  = - u[:, rank_in_data]
    vh[rank_in_data, :] = -vh[rank_in_data, :]

    return None

job_name = "xfam"

# [[[ LOAD U DATA ]]]
# Load upstream data...
u     = np.load("u.seq.npy")
s     = np.load("s.seq.npy")
vh    = np.load("vh.seq.npy")
len_res = np.load("len_res.npy")
len_seq = np.load("len_seq.npy")

# Allow positive value that fits the distance in nature (optinal)...
reverse_sign(u, vh, 1, index_from_zero = False)
rev_list = [7]
for rev_i in rev_list: reverse_sign(u, vh, rev_i, index_from_zero = False)

# Calculate the coefficients...
c = np.matmul(np.diag(s), vh)

# Standardize u and c and assign units...
u_ave = np.sqrt(u.shape[0])
c = c / u_ave

# Define a series of rotation...
rotations = [
    [2, 3, 30.9-4.8],
    [2, 4, -8.1],
    [3, 4, 38.3],
    [3, 5, 14.8],
    [3, 6, -10.2],
    [4, 5, -24.4],
    [4, 6, 3.3],
    [4, 11, -41.8],
    [6,  8, -41.0],
    [6,  7,  25.1],
    [8,  2, -10.2],
    [8,  7, 18.7],
    [8,  9,-11.1],
    [9,  6,-33.5],
    [9,  7,-14.3],
    [9, 10,-22.9],
    [7, 11,-18.1],
    [7,  5,-22.8],
    [7,  4,-11.4],
    [2,  7,  5.6],
]
for i, (x, y, _) in enumerate(rotations):
    if x in rev_list: rotations[i][2] *= -1
    if y in rev_list: rotations[i][2] *= -1
disp_index = -1    # 0-based Python convention
if len(rotations): rank1_last, rank2_last = rotations[disp_index][0:0 + 2]
for rank1, rank2, theta in rotations:
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)


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


# [[[ CONVERT TO GRADIENT ]]]
# [IMPROVE] Wrap this?
#  0. gray=0.0000, (r,g,b)=(0.5020,0.0000,0.0000), #800000 = 128   0   0
#  1. gray=0.2500, (r,g,b)=(1.0000,0.0000,0.0000), #ff0000 = 255   0   0
#  2. gray=0.4500, (r,g,b)=(1.0000,1.0000,1.0000), #ffffff = 255 255 255
#  3. gray=0.5000, (r,g,b)=(0.7569,1.0000,0.7569), #c1ffc1 = 193 255 193
#  4. gray=0.5500, (r,g,b)=(1.0000,1.0000,1.0000), #ffffff = 255 255 255
#  5. gray=0.7500, (r,g,b)=(0.0000,0.0000,1.0000), #0000ff =   0   0 255
#  6. gray=1.0000, (r,g,b)=(0.0000,0.0000,0.5020), #000080 =   0   0 128
resol  = 10
c_list = []
c_list.extend( cs.linear_gradient("800000", "ff0000", n = (-5 - (-10)) * resol ) )
c_list.extend( cs.linear_gradient("ff0000", "ffffff", n = (-1 - (-5) ) * resol ) )
c_list.extend( cs.linear_gradient("ffffff", "c1ffc1", n = (0  - (-1) ) * resol ) )
c_list.extend( cs.linear_gradient("c1ffc1", "ffffff", n = (1  - ( 0) ) * resol ) )
c_list.extend( cs.linear_gradient("ffffff", "0000ff", n = (5  - ( 1) ) * resol ) )
c_list.extend( cs.linear_gradient("0000ff", "000080", n = (10 - ( 5) ) * resol ) )

def guess_color(c_list, r):
    ''' r is ratio
    '''
    l = len(c_list)
    ## print(r, l, l * r)
    return c_list[int(l * r) - 1]


# [[[ MAP SEQI to COLOR ]]]
# [IMPROVE] Wrap this?
rank = 7
dmat = pr.utils.array2tril(u[:, rank - 1], len_seq, offset = -1)
dmat_full = dmat + dmat.T
column_mean_dmat = np.nanmean(dmat_full, axis = 0, keepdims = False)
mean_min = np.min(column_mean_dmat)
mean_max = np.max(column_mean_dmat)
mean_abs = np.max([np.abs(mean_min), np.abs(mean_max)])

# Ratios
rats = (column_mean_dmat - (-mean_abs)) / (2 * mean_abs)
seqi_to_color_dict = { nseqi + i//len_res: guess_color(c_list, rat) for i, rat in enumerate(rats.tolist()) }


# [[[ LOAD DATABASE ]]]
# Specify the rigid framework...
fl_fwk = 'fwk.dat'
fwk    = pr.utils.read_file(fl_fwk)
fwk    = [ [int(i) // len_res, int(j) // len_res] for i, j in fwk ]

# Specify the backbone atoms to select...
backbone_select = "name " + "+".join( backbone )




# ======== MOLECLUE REFERENCE =========
# Pick a PDB...
pdb, chain = "6PS6", "A"
drc = "pdb"
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
        seqi_fwk = nseqi + i
        resi_fwk = seqi_to_resi_dict[seqi_fwk]
        if resi_fwk is None: continue
        fwk_resi.append(resi_fwk)
fwk_select = {}
fwk_select[f"{pdb}_{chain}"] = '+'.join( [ str(i) for i in fwk_resi ] )

# Map resi to color...
resi_to_color_dict = {}
resi_to_color_dict = { seqi_to_resi_dict[seqi] : color for seqi, color in seqi_to_color_dict.items() }


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

for resi, color in resi_to_color_dict.items():
    pm(f"set cartoon_color, 0x{color}, %{entry} and resi {resi}")




# ======== MOLECLUE MOBILE =========
pdb, chain = "4LDE", "A"
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
        seqi_fwk = nseqi + i
        resi_fwk = seqi_to_resi_dict[seqi_fwk]
        if resi_fwk is None: continue
        fwk_resi.append(resi_fwk)

fwk_select[f"{pdb}_{chain}"] = '+'.join( [ str(i) for i in fwk_resi ] )

# Map resi to color...
resi_to_color_dict = {}
resi_to_color_dict = { seqi_to_resi_dict[seqi] : color for seqi, color in seqi_to_color_dict.items() }

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

for resi, color in resi_to_color_dict.items():
    pm(f"set cartoon_color, 0x{color}, %{entry} and resi {resi}")

# Align...
pm(f"super {mobile}, {target}")


pm("orient")
pm("window size, 1500, 1500")

input()
