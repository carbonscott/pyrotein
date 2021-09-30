#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pyrotein as pr
import givens as gv
from display import plot_dmat, plot_singular, plot_left_singular, \
                    plot_coeff, select_items, plot_blankcoeff
import multiprocessing as mp
from loaddata import load_gpcrdb_xlsx, load_md_xlsx
import colorsimple as cs
from itertools import product, permutations
import os


# [[[ Experimental ]]]
def reverse_sign(u, vh, rank, index_from_zero = True):
    # Comply with the convention (1-based index)
    rank_in_data = rank if index_from_zero else rank - 1

    # Reverse sign...
    u[:, rank_in_data]  = - u[:, rank_in_data]
    vh[rank_in_data, :] = -vh[rank_in_data, :]

    return None


# Set job name...
job_name = "xfam-loop"
aim      = "class+md"
fl_path  = f"{job_name}.{aim}.c"

# Specify chains to process...
fl_chain = "gpcrdb.all.xlsx"
sheet    = f"xfam"
lines    = load_gpcrdb_xlsx(fl_chain, sheet = sheet, splitchain = True)

# Specify md trajectories to process...
sheet    = f"md"
md_lines = load_md_xlsx(fl_chain, sheet = sheet)

# Load upstream data...
b, e  = 0, 20 + 1
u       = np.load(f"{job_name}.u.seq.trunc.{b:02d}to{e:02d}.npy")
s       = np.load(f"{job_name}.s.seq.trunc.{b:02d}to{e:02d}.npy")
vh      = np.load(f"{job_name}.vh.seq.trunc.{b:02d}to{e:02d}.npy")
len_res = np.load(f"{job_name}.len_res.npy")
len_seq = np.load(f"{job_name}.len_seq.npy")

# Allow positive value that fits the distance in nature (optinal)...
reverse_sign(u, vh, 1, index_from_zero = False)
rev_list = [3, 4]
for rev_i in rev_list: reverse_sign(u, vh, rev_i, index_from_zero = False)

# Calculate the weighted coefficients...
# The following is a computational equivalence to c = RMS(u @ s, axis = 0) @ vh
c = np.matmul(np.diag(s), vh)
u_rms = 1.0 / np.sqrt(u.shape[0])
c     = c * u_rms

# Define a series of rotation...
rotations = [
    [  2,  3, -48.8],
    [  2,  5, -14.5],
    [  3,  5, -43.4],
    [  5,  7,  43.4],
    [  8,  9,  43.4],
    [  3,  4,  19.2],
    [  2,  4,  13.4],
]
for i, (x, y, _) in enumerate(rotations):
    if x in rev_list: rotations[i][2] *= -1
    if y in rev_list: rotations[i][2] *= -1
disp_index = -1    # 0-based Python convention
if len(rotations):
    rank1_last, rank2_last = rotations[disp_index][0:0 + 2]
for rank1, rank2, theta in rotations:
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)

## rank1_last, rank2_last = 3, 5




# [[[ MD Find the tar_seq for structural alignment ]]]
# Read the MD trajectory recrods...
drc_md = "mdtraj"
md_recs = []
for i, line in enumerate(md_lines):
    pdb, chain, frame_b, frame_e, frame_s, chain_orig = line[1:7]
    for i in range(frame_b, frame_e, frame_s):
        md_rec = (pdb, chain, f"{i:02d}", chain_orig)
        md_recs.append(md_rec)

# Read the sequence alignment result...
fl_aln   = f"xfam.step4.psa.fil.fasta"
seq_dict = pr.fasta.read(fl_aln)

# Obtain the consensus sequence (super seq)...
tally_dict = pr.fasta.tally_resn_in_seqs(seq_dict)
super_seq  = pr.fasta.infer_super_seq(tally_dict)

# The default ones are not good enough
nseqi = 637
cseqi = 1371

# Define atoms used for distance matrix analysis...
backbone = ["N", "CA", "C", "O"]
len_res = len(backbone)
len_seq = (cseqi - nseqi + 1) * len_res

# Accumulate distance matices as lower triangluar matrix...
len_lower_tri = (len_seq * len_seq - len_seq) // 2
dmats = np.zeros((len(md_recs), len_lower_tri))

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
for _, (b, e) in labels.items(): mask_loop[b:e, :] = 1.0

# Download the entrywise mean for dmats...
dmats_entrywise_mean = np.load(f"{job_name}.dmats.entrywise_mean.npy").reshape(-1)

# Go through each record...
for i, md_rec in enumerate(md_recs):
    pdb, chain, frame, chain_orig = md_rec

    fl_pdb    = f"{pdb}.{frame}.mod.pdb"
    pdb_path  = os.path.join(drc_md, fl_pdb)
    atoms_pdb = pr.atom.read(pdb_path)

    # Create a lookup table for this pdb...
    atom_dict = pr.atom.create_lookup_table(atoms_pdb)

    # Obtain the chain to process...
    chain_dict = atom_dict[chain]

    # Obtain seq string for the current chain...
    tar_seq = seq_dict[f"{pdb.lower()}_{chain_orig}"]

    # Obtain xyzs
    entry = f"{pdb}_{chain}"
    print(entry)
    xyzs = pr.atom.extract_xyz_by_seq(backbone, chain_dict, tar_seq, nseqi, cseqi)
    xyzs[mask_loop == 0] = np.nan
    ## xyzs[mask_loop == 0] = 0.0

    # Calculate distance matrix...
    dmat = pr.distance.calc_dmat(xyzs, xyzs)
    dmat1d = pr.utils.mat2tril(dmat, offset = -1)

    # Patch nan based on column mean obtained prior to SVD...
    nan_pos = np.isnan(dmat1d)
    dmat1d[nan_pos] = dmats_entrywise_mean[nan_pos]
    dmats[i, :] = dmat1d


c_md = u.T @ dmats.T
c_md = c_md * u_rms

c_composite = np.hstack((c, c_md))

# [[[ Plot style ]]]
# Write this get_pointstyle function for different coloring purposes
# The point is that each dot should be customizable
# Delineate the input (degree of freedom) and output
# f(cls, act, res) in this case that defines the look of a point
# The following function illustrates the mapping, however, global variable 
# can resolve some memory overhead.  So in practice, it's not wrapped in a 
# function.  
#
# ```
# def get_pointstyle(cls, act, res):
#     ''' A style is defined by three variables: class, activation, resolution.
#     '''
#     # Color is defined by class...
#     clr_dict = {
#         "A (Rhodopsin)" : "#FFB8B8" ,   # Low key red
#         "B1 (Secretin)" : "#00A600" ,
#         "C (Glutamate)" : "#0080FF" ,
#         "B2 (Adhesion)" : "gray"    ,
#         "F (Frizzled)"  : "#AB00FF" ,
#     }
# 
#     # Shape is defined by activation...
#     shape_dict = {
#         "Intermediate" : 8,
#         "Inactive"     : 6,
#         "Active"       : 2,
#     }
# 
#     # Pointsize is defined by resolution...
#     pointsize_dict = {
#         "low" : 0.3,
#         "hi"  : 1.0,
#     }
# 
#     # Fetch values from the rules based on the input...
#     clr       = clr_dict      [cls]
#     shape     = shape_dict    [act]
#     pointsize = pointsize_dict[res]
# 
#     cmd = f"u 1:2 w p pt {shape} pointsize {pointsize} lw 1.0 lc rgb '{clr}' notitle"
# 
#     return cmd
# ```

# PART 1: Experimental data
# Associate STYLE ELEMENT to ENTRY CHARACTERISTICS
# Color is defined by class
clr_dict = {
    "A (Rhodopsin)" : "#FFB8B8" ,   # Low key red
    "B1 (Secretin)" : "#00A600" ,
    "C (Glutamate)" : "#0080FF" ,
    ## "B2 (Adhesion)" : "gray"    ,
    "F (Frizzled)"  : "#AB00FF" ,
}

# Shape is defined by activation
shape_dict = {
    "Intermediate" : 8,
    "Inactive"     : 6,
    "Active"       : 2,
}

# Pointsize is defined by resolution
pointsize_dict = {
    "low" : 0.3,
    "hi"  : 1.0,
}
# Resolution cutoff...
res_cutoff = 3.5

# Compose the plotting command based on style element
def get_pointstyle(clr, shape, pointsize):
    ''' A style is defined by three variables: class, activation, resolution.
    '''
    return f"u 1:2 w p pt {shape} pointsize {pointsize} lw 1.0 lc rgb '{clr}' notitle"

# Plot style for each input...
gps_dict = {}
for i, line in enumerate(lines[:]):
    pdb, chain = line[7], line[10]

    # Fetch characteristics from each entry...
    cls = line[4]
    act = line[11]
    res = float(line[9])
    res_label = "hi" if res < res_cutoff else "low"

    # Fetch style element from the rules based for this input...
    # You don't have to plot every record in that column
    if not cls in clr_dict            : continue
    if not act in shape_dict          : continue
    if not res_label in pointsize_dict: continue
    clr       = clr_dict      [cls]
    shape     = shape_dict    [act]
    pointsize = pointsize_dict[res_label]

    # Define the style for this input...
    entry_name = f"{pdb}_{chain}"
    gps_dict[entry_name] = { "style" : get_pointstyle(clr, shape, pointsize), 
                             "entry" : [i] }


# PART 2: MD data
# Continue to use gps_dict
md_clr_dict = {
    "1u19" : "blue",
    "3pqr" : "red",
}

pdb, chain, frame, chain_orig = md_recs[0]
gps_dict[f"md.{pdb}"] = {
    "style" : f"u 1:2 w l lc rgb '{md_clr_dict[pdb]}' notitle", 
    "entry" : [i + len(lines) for i, _ in enumerate(md_recs)]
}
gps_dict[f"md.{pdb}.head"] = {
    "style" : f"u 1:2 w p pt 6 ps 0.5 lc rgb '{md_clr_dict[pdb]}' notitle", 
    "entry" : [len(lines)],
}
gps_dict[f"md.{pdb}.tail"] = {
    "style" : f"u 1:2 w p pt 7 ps 0.5 lc rgb '{md_clr_dict[pdb]}' notitle", 
    "entry" : [len(lines) + len(md_recs)-1],
}


cs.color_table(gps_dict)

# [[[ PLOT ]]]
# Create labels...
entry_dict = { f"{v[7]}_{v[10]}" : i for i, v in enumerate(lines) }

# Selectively plot entries...
pdb_list = [ 
           ]

entry_fil_dict = { k : entry_dict[k] for k in pdb_list if k in entry_dict }
## entry_fil_dict = entry_dict.copy()

cmds = [
           f"set xzeroaxis",
           f"set yzeroaxis",
       ]

top = 20 + 1


def run(rank1, rank2):
    offset = "1.8,0.0"
    cmds.append(f"set size ratio -1")
    return plot_coeff(c_composite, rank1, rank2,
                      lbl = gps_dict,
                      ## label = True, 
                      label = False, 
                      label_dict = entry_fil_dict,
                      lbl_fontsize = 6,

                      ## xrange = (-0.30,-0.20),
                      ## yrange = (-0.20,-0.10),

                      offset = offset, 
                      rot = 0,
                      height = 2.873,
                      width = 2.873,
                      fontsize = 12,
                      pointsize = 1.0,
                      linewidth = 1.0,
                      fl_path = fl_path,
                      fl_postfix = f'',
                      is_rug = False,
                      index_from_zero = False,
                      cmds = cmds)


if 1:
    if 0:
        # Visualize singular values...
        top = 10
        plot_singular(s, top = top, log = True, 

                      width = 6,
                      height = 5,
                      fontsize = 32,

                      ## width = 4,
                      ## height = 3,
                      ## fontsize = 16,

                      linewidth = 1.0,
                      ticscale  = 2.0,
                      fl_path = fl_path,
                      index_from_zero = False)
    if 0:
        rank1, rank2 = 2, 3
        run(rank1, rank2)

    if 0:
        rank1_adjust, rank2_adjust = rank1_last, rank2_last
        num_job = 10
        start_index = 1
        top = 20 + 1
        job_rng = range(1, top)
        ## job_ids = product(job_rng, job_rng)
        job_ids = filter(lambda x: x[0] < x[1], product([rank1_adjust], job_rng))
        with mp.Pool(num_job) as proc:
            proc.starmap( run, job_ids )
        job_ids = filter(lambda x: x[0] < x[1], product([rank2_adjust], job_rng))
        with mp.Pool(num_job) as proc:
            proc.starmap( run, job_ids )
        job_ids = filter(lambda x: x[0] < x[1], product(job_rng, [rank1_adjust]))
        with mp.Pool(num_job) as proc:
            proc.starmap( run, job_ids )
        job_ids = filter(lambda x: x[0] < x[1], product(job_rng, [rank2_adjust]))
        with mp.Pool(num_job) as proc:
            proc.starmap( run, job_ids )

    if 1:
        start_index = 1
        top = 10 + 1
        job_rng = range(1, top)
        job_ids = filter(lambda x: x[0] < x[1], permutations(job_rng, 2))
        num_job = 4
        if __name__ == "__main__":
            with mp.Pool(num_job) as proc:
                proc.starmap( run, job_ids )



    if 0:
        offset = "1.8,0.0"
        start_index = 1
        for j in range(start_index, top):
            ## for i in range(start_index + 1, top):
            for i in range(1, top):
                rank1, rank2 = j, i
                plot_coeff(c, rank1, rank2,
                           lbl = gps_dict,
                           ## label = True, 
                           label = False, 
                           label_dict = entry_dict,
                           lbl_fontsize = 10,
                           ## xrange = (-2.0, 2.0),
                           ## yrange = (-2.0, 2.0),
                           offset = offset, 
                           rot = 0,
                           height = 6,
                           width = 6,
                           fontsize = 25,
                           pointsize = 2.0,
                           linewidth = 1.0,
                           fl_path = fl_path,
                           fl_postfix = f'',
                           index_from_zero = False,
                           cmds = cmds)

    if 0:
        cmds.append(f"set xtics 0.5")
        offset = "1.8,0.0"
        start_index = 3
        for j in range(start_index, top):
            rank1, rank2 = j, 1
            plot_coeff(c, rank1, rank2,
                       lbl = gps_dict,
                       label = False, 
                       label_dict = entry_dict,
                       lbl_fontsize = 10,
                       xrange = (-2.5, 2.5),
                       yrange = (23.6, 22.6),
                       offset = offset, 
                       rot = 0,
                       height = 6,
                       width = 6,
                       fontsize = 25,
                       pointsize = 2.0,
                       linewidth = 1.0,
                       fl_path = fl_path,
                       fl_postfix = f'',
                       index_from_zero = False,
                       cmds = cmds)

    if 0:
        offset = "1.8,0.0"
        cmds.append("set size ratio -1")
        rank1, rank2 = 4, 2
        plot_coeff(c, rank1, rank2,
                   lbl = gps_dict,
                   label = True, 
                   ## label_dict = entry_dict,
                   label_dict = entry_fil_dict,
                   lbl_fontsize = 10,
                   xrange = (-2.0, 2.0),
                   yrange = (-2.0, 2.0),
                   offset = offset, 
                   rot = 0,
                   height = 6,
                   width = 6,
                   fontsize = 25,
                   pointsize = 2.0,
                   linewidth = 1.0,
                   fl_path = fl_path,
                   fl_postfix = f'',
                   index_from_zero = False,
                   cmds = cmds)

    if 0:
        offset = "1.8,0.0"
        cmds.append("set encoding utf8")
        ## cmds.append("set key top left")
        start_index = 1
        rank1_adjust, rank2_adjust = rank1_last, rank2_last
        for i in [rank1_adjust, rank2_adjust]:
            for j in range(start_index, top):
                rank1, rank2 = i, j
                plot_coeff(c, rank1, rank2,
                           lbl = gps_dict,
                           ## label = True, 
                           label = False, 
                           label_dict = entry_fil_dict,
                           lbl_fontsize = 10,
                           xrange = (-3.0, 3.0),
                           yrange = (-3.0, 3.0),
                           offset = offset, 
                           rot = 0,
                           height = 6,
                           width = 6,
                           fontsize = 25,
                           pointsize = 2.0,
                           linewidth = 1.0,
                           fl_path = fl_path,
                           fl_postfix = f'',
                           index_from_zero = False,
                           cmds = cmds)

                rank1, rank2 = j, i
                plot_coeff(c, rank1, rank2,
                           lbl = gps_dict,
                           ## label = True, 
                           label = False, 
                           label_dict = entry_fil_dict,
                           lbl_fontsize = 10,
                           xrange = (-3.0, 3.0),
                           yrange = (-3.0, 3.0),
                           offset = offset, 
                           rot = 0,
                           height = 6,
                           width = 6,
                           fontsize = 25,
                           pointsize = 2.0,
                           linewidth = 1.0,
                           fl_path = fl_path,
                           fl_postfix = f'',
                           index_from_zero = False,
                           cmds = cmds)
