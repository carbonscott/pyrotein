#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pyrotein as pr
from display import plot_rmsd_dmat, calc_rigid_framework
import pprint

# Set job name...
job_name = "xfam-loop"

# load the data...
rmsd_dmat = np.load(f"{job_name}.rmsd_dmat.seq.npy")
len_res   = np.load(f"{job_name}.rmsd_len.seq.npy")
nseqi     = np.load(f"{job_name}.nseqi.seq.npy")
cseqi     = np.load(f"{job_name}.cseqi.seq.npy")
num_seq   = cseqi - nseqi + 1

fl_rmsd_dmat  = f"{job_name}.rmsd"

# Create labels...
# +1 to make it a right-end closed
labels = {'H8': [722, 735],
          'TM1': [0, 34],
          'TM2': [63, 96],
          'TM3': [158, 194],
          'TM4': [227, 254],
          'TM5': [316, 357],
          'TM6': [591, 631],
          'TM7': [680, 708]}
# Increment the right index so that it is right-included
for k, v in labels.items(): labels[k][1] += 1

# Leave rooms for all atoms in a residue
for k, v in labels.items(): labels[k] = [ (i) * len_res for i in v ]

# Read the initial framework...
# +1 to right-include the index has been considered in the fl_fwk
fl_fwk = f"fwk.{job_name}.dat"
seqi_fwk_list = pr.utils.read_file(fl_fwk, str_to_num = True, num_type = int)
fwk_list = calc_rigid_framework(rmsd_dmat, seqi_fwk_list, num_seq, len_res, min_size = 5, min_mean = 0.2, epsilon = 0.0001)


# Define a colorscheme...
# Colorscheme is inspired by from this publication (DOI: 10.1093/nar/gkw555) from Zhong Ren
pal = "set palette defined ( 0 'seagreen', 0.1 'white', 0.5 'blue', 1 'navy' )"
plot_rmsd_dmat( rmsd_dmat, 
                fl_rmsd_dmat, 
                fwk_list = fwk_list,
                fwk_linewidth = 3,
                curve_linewidth = 3,
                lbl = labels,
                width = 10,
                height = 12,
                fontsize = 29,
                lbl_fontsize = 29,
                linewidth = 1.0,
                palette = pal,
                intst_max = "*",
              )
