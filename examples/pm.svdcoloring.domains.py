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

comb = {
    "c2 A" : { "pdb1"             : "1u19_A"      ,
               "pdb2"             : "3cap_A"      ,
               "rank"             : 2             ,
               "domain_div"       : "TM3"         ,
               "domain_intervals" : 2             ,
               "domain_hl"        : 2             ,
               "drc_pdb"          : "pdb.numbered",
               "postfix"          : "_GPCRDB"     , },


    "c3 B1" : { "pdb1"             : "5VEW_A"      ,
                "pdb2"             : "6X18_R"      ,
                "rank"             : 3             ,
                "domain_div"       : "TM6"         ,
                "domain_intervals" : 3             ,
                "domain_hl"        : 1             ,
                "drc_pdb"          : "pdb",
                "postfix"          : ""     , },
    "c2 B1" : { "pdb1"             : "5VEW_A"      ,
                "pdb2"             : "6X18_R"      ,
                "rank"             : 2             ,
                "domain_div"       : "TM6"         ,
                "domain_intervals" : 3             ,
                "domain_hl"        : 1             ,
                "drc_pdb"          : "pdb",
                "postfix"          : ""     , },


    "c2 F" : { "pdb1"             : "4QIM_A"      ,
               "pdb2"             : "6XBK_R"      ,
               "rank"             : 2             ,
               "domain_div"       : "TM6"         ,
               "domain_intervals" : 3             ,
               "domain_hl"        : 1             ,
               "drc_pdb"          : "pdb",
               "postfix"          : ""     , },
    "c3 F" : { "pdb1"             : "4QIM_A"      ,
               "pdb2"             : "6XBK_R"      ,
               "rank"             : 3             ,
               "domain_div"       : "TM6"         ,
               "domain_intervals" : 3             ,
               "domain_hl"        : 1             ,
               "drc_pdb"          : "pdb",
               "postfix"          : ""     , },
}

choice = "c2 A"



# [[[ Hide from users ]]]

# Fetch values
pdb1 = comb[choice]["pdb1"].split("_")
pdb2 = comb[choice]["pdb2"].split("_")

rank = comb[choice]["rank"]
domain_div = comb[choice]["domain_div"]
domain_intervals = comb[choice]["domain_intervals"] # Divide domain into chunks
domain_hl = comb[choice]["domain_hl"]  # Highlight one chunk

drc_pdb = comb[choice]["drc_pdb"]
postfix = comb[choice]["postfix"]

# Format pdb
pdb1[0] = pdb1[0].lower()
pdb2[0] = pdb2[0].lower()

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

# Derive the domains to divide...
domain_b, domain_e = labels[domain_div]
domain_len = domain_e - domain_b
domain_interval_len = domain_len // domain_intervals
chunk_b = domain_b + (domain_hl - 1) * domain_interval_len    # index from 1
chunk_e = domain_b + domain_hl * domain_interval_len
box_dict = {
    "chunk_ref" : [chunk_b, chunk_e],
}


# Create mask...
# +1 offset is to deal with numpy indexing which is not right-inclusive 
labels_aux = {}
for k, v in labels.items(): labels_aux[k] = [ (i) * len_res for i in v ]
for _, (b, e) in labels_aux.items(): mask_loop[b:e+1, :] = 1

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


# [[[ LOAD U DATA ]]]
# Load upstream data...
b = 0
e = 20 + 1
u     = np.load(f"{job_name}.u.seq.trunc.{b:02d}to{e:02d}.npy")
s     = np.load(f"{job_name}.s.seq.trunc.{b:02d}to{e:02d}.npy")
vh    = np.load(f"{job_name}.vh.seq.trunc.{b:02d}to{e:02d}.npy")
len_res = np.load(f"{job_name}.len_res.npy")
len_seq = np.load(f"{job_name}.len_seq.npy")

def reverse_sign(u, vh, rank, index_from_zero = True):
    # Comply with the convention (1-based index)
    rank_in_data = rank if index_from_zero else rank - 1

    # Reverse sign...
    u[:, rank_in_data]  = - u[:, rank_in_data]
    vh[rank_in_data, :] = -vh[rank_in_data, :]

    return None

# Allow positive value that fits the distance in nature (optinal)...
reverse_sign(u, vh, 1, index_from_zero = False)
rev_list = [3, 4]
for rev_i in rev_list: reverse_sign(u, vh, rev_i, index_from_zero = False)

# Calculate the coefficients...
c = np.matmul(np.diag(s), vh)

# Standardize u and c and assign units...
u_ave = np.sqrt(u.shape[0])
c = c / u_ave

# Define a series of rotation...
# NEED THE LATEST ONE!!!
rotations = [
    [  2,  3, -48.8],
    [  2,  5, -14.5],
    [  3,  5, -43.4],
    [  5,  7,  43.4],
    [  8,  9,  43.4],
    [  3,  4,  19.2],
    [  2,  4,  13.4],
    [  5,  6, -21.3],
]
for i, (x, y, _) in enumerate(rotations):
    if x in rev_list: rotations[i][2] *= -1
    if y in rev_list: rotations[i][2] *= -1
disp_index = -1    # 0-based Python convention
if len(rotations): rank1_last, rank2_last = rotations[disp_index][0:0 + 2]
for rank1, rank2, theta in rotations:
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)


# [[[ Pairwise Product ]]]
rank_in_data = rank - 1
u_rank = u[:, rank_in_data]

# Find the full range...
frac = 1.0
bound = np.nanmax(np.abs([np.nanmin(u_rank), np.nanmax(u_rank)]))
intst_min = -bound * frac
intst_max =  bound * frac

# Convert `u` at position `rank` into a lower triangular matrix...
u_rank = pr.utils.array2tril(u_rank, len_seq, offset = -1)

# Restore dmat to a full matrix for visualization...
u_rank = u_rank + u_rank.T


# [[[ CONVERT TO GRADIENT ]]]
# [IMPROVE] Wrap this?
#  0. gray=0.0000, (r,g,b)=(0.5020,0.0000,0.0000), #800000 = 128   0   0
#  1. gray=0.2500, (r,g,b)=(1.0000,0.0000,0.0000), #ff0000 = 255   0   0
#  2. gray=0.4500, (r,g,b)=(1.0000,1.0000,1.0000), #ffffff = 255 255 255
#  3. gray=0.5000, (r,g,b)=(0.7569,1.0000,0.7569), #c1ffc1 = 193 255 193
#  4. gray=0.5500, (r,g,b)=(1.0000,1.0000,1.0000), #ffffff = 255 255 255
#  5. gray=0.7500, (r,g,b)=(0.0000,0.0000,1.0000), #0000ff =   0   0 255
#  6. gray=1.0000, (r,g,b)=(0.0000,0.0000,0.5020), #000080 =   0   0 128

resol  = 100
c_list = []
c_list.extend( cs.linear_gradient("800000", "ff0000", n = (-5 - (-10)) * resol, sym = "0x") )
c_list.extend( cs.linear_gradient("ff0000", "ffffff", n = (-1 - (-5) ) * resol, sym = "0x") )
c_list.extend( cs.linear_gradient("ffffff", "c1ffc1", n = (0  - (-1) ) * resol, sym = "0x") )
c_list.extend( cs.linear_gradient("c1ffc1", "ffffff", n = (1  - ( 0) ) * resol, sym = "0x") )
c_list.extend( cs.linear_gradient("ffffff", "0000ff", n = (5  - ( 1) ) * resol, sym = "0x") )
c_list.extend( cs.linear_gradient("0000ff", "000080", n = (10 - ( 5) ) * resol, sym = "0x") )

def guess_color(c_list, v, vmin = -0.01, vmax = 0.01):
    ''' r is ratio
    '''
    l = len(c_list)
    ## print(r, l, l * r)

    # Handle color out of the range...
    if v > vmax: return c_list[-1]
    if v < vmin: return c_list[0]

    # Handle color in range...
    r = (v - vmin) / (vmax - vmin)
    ## print(r, v)

    return c_list[int((l-1) * r)]


# [[[ BOX RANGE ]]]
box_atom_b, box_atom_e = [ i * len_res for i in box_dict["chunk_ref"] ]


# [[[ MAP SEQI to COLOR ]]]
# [IMPROVE] Wrap this?
column_mean_dmat = np.nanmean(u_rank[box_atom_b:box_atom_e,:], axis = 0, keepdims = False)

seqi_to_color_dict = {}
for i, seqi in enumerate(range(nseqi, cseqi + 1)):
    val_ave = np.mean(column_mean_dmat[i * len_res : (i+1) * len_res])
    seqi_to_color_dict[seqi] = guess_color(c_list, val_ave, vmin = -0.005, vmax = 0.005)



# [[[ LOAD DATABASE ]]]
# Specify the backbone atoms to select...
backbone_select = "name " + "+".join( backbone )


# ======== MOLECLUE REFERENCE =========
# Pick a PDB...
pdb, chain = pdb1
entry = f"{pdb}{postfix}"
fl_pdb     = f"{entry}.pdb"
pdb_path   = os.path.join(drc_pdb, fl_pdb)
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

# Map resi to color...
# WRONG
resi_to_color_dict = {}
for seqi, color in seqi_to_color_dict.items():
    resi = seqi_to_resi_dict[seqi]

    # Skip those gaps...
    if resi == None: continue

    # Assign color...
    resi_to_color_dict[resi] = color


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
pdb_path   = os.path.join(drc_pdb, f"{entry}.pdb")
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

seqi_b, seqi_e = box_dict["chunk_ref"]
resi_b, resi_e = seqi_to_resi_dict[seqi_b + nseqi], seqi_to_resi_dict[seqi_e + nseqi]

# Auto bypass resi that is None
i = 1
while resi_b is None and seqi_b + nseqi + i <= cseqi: 
    resi_b = seqi_to_resi_dict[seqi_b + nseqi + i]
    i += 1
i = 1
while resi_e is None and seqi_e + nseqi - i >= nseqi: 
    resi_e = seqi_to_resi_dict[seqi_e + nseqi - i]
    i += 1

pm(f"set cartoon_color, yellow, %{entry} and resi {resi_b}-{resi_e}")



# ======== MOLECLUE MOBILE 1 =========
pdb, chain = pdb2
entry      = f"{pdb}{postfix}"
pdb_path   = os.path.join(drc_pdb, f"{entry}.pdb")
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

seqi_b, seqi_e = box_dict["chunk_ref"]
resi_b, resi_e = seqi_to_resi_dict[seqi_b + nseqi], seqi_to_resi_dict[seqi_e + nseqi]

# Auto bypass resi that is None
i = 1
while resi_b is None and seqi_b + nseqi + i <= cseqi: 
    resi_b = seqi_to_resi_dict[seqi_b + nseqi + i]
    i += 1
i = 1
while resi_e is None and seqi_e + nseqi - i >= nseqi: 
    resi_e = seqi_to_resi_dict[seqi_e + nseqi - i]
    i += 1

pm(f"set cartoon_color, yellow, %{entry} and resi {resi_b}-{resi_e}")

# Align...
pm(f"super {mobile}, {target}")


pm("orient")
pm("window size, 1500, 1500")

input()
