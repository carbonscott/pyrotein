#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pyrotein as pr
from loaddata import load_xlsx
from display import plot_dmat
import pprint


# load the data...
rmsd_dmat = np.load("rmsd_dmat.seq.npy")
len_res   = np.load("rmsd_len.seq.npy")

fl_rmsd_dmat  = "rmsd"

# Define a colorscheme...
# Colorscheme is inspired by from this publication (DOI: 10.1093/nar/gkw555) from Zhong Ren
pal = "set palette defined ( 0 'seagreen', 0.1 'white', 0.5 'blue', 1 'navy' )"

# Create labels...
# +1 to make it a right-end closed
## labels = {
##     "TM1" : [ 36,  60 + 1],
##     "TM2" : [ 69,  97 + 1],
##     "TM3" : [104, 136 + 1],
##     "TM4" : [146, 171 + 1],
##     "TM5" : [189, 221 + 1],
##     "TM6" : [369, 397 + 1],
##     "TM7" : [406, 427 + 1],
##     "H8"  : [432, 441 + 1],
## }
## for k, v in labels.items(): labels[k] = [ (i - nterm) * len_res for i in v ]

plot_dmat( rmsd_dmat, 
           fl_rmsd_dmat, 
           ## lbl = labels,
           width = 10,
           height = 12,
           fontsize = 29,
           lbl_fontsize = 29,
           linewidth = 1.0,
           palette = pal,
           intst_max = "*",
         )


