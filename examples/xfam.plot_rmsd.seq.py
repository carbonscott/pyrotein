#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pyrotein as pr
from display import plot_rmsd_dmat
import pprint

# Set job name...
job_name = "xfam"

# load the data...
rmsd_dmat = np.load("rmsd_dmat.seq.npy")
len_res   = np.load("rmsd_len.seq.npy")
nseqi     = np.load("nseqi.seq.npy")
cseqi     = np.load("cseqi.seq.npy")

fl_rmsd_dmat  = f"{job_name}.rmsd"

# Define a colorscheme...
# Colorscheme is inspired by from this publication (DOI: 10.1093/nar/gkw555) from Zhong Ren
pal = "set palette defined ( 0 'seagreen', 0.1 'white', 0.5 'blue', 1 'navy' )"

# Create labels...
# +1 to make it a right-end closed
labels = {'H8': [620, 633],
          'TM1': [0, 33],
          'TM2': [47, 79],
          'TM3': [115, 151],
          'TM4': [182, 209],
          'TM5': [277, 317],
          'TM6': [503, 542],
          'TM7': [584, 612]}
for k, v in labels.items(): labels[k] = [ (i) * len_res for i in v ]

plot_rmsd_dmat( rmsd_dmat, 
                fl_rmsd_dmat, 
                pop_bin_cap = 50,
                fwk_mid = 2.072870590257889,
                fwk_tol = 0.5,
                fwk_minsize = 10,
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

