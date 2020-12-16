#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pyrotein as pr
from display import plot_dmat, plot_singular, plot_left_singular, plot_coeff

# Specify the range of atoms from rhodopsin...
nterm = 1
cterm = 348
backbone = ["N", "CA", "C", "O"]
length_backbone = (cterm - nterm + 1) * len(backbone)

# Load upstream data...
lines = np.load("lines.npy").tolist()
dmats = np.load("dmats.npy")
u     = np.load("u.npy")
s     = np.load("s.npy")
vh    = np.load("vh.npy")
c     = np.load("c.npy")

if 0:
    # Visualize singular values...
    plot_singular(s, log = True)

if 0:
    # Visualize a u matrix...
    for rank in range(1,10):
        plot_left_singular(u, rank, 
                              length_mat = length_backbone, 
                              frac = 1.0,
                              binning = 2)

if 1:
    # Create labels...
    labels = ['_'.join(i) for i in lines]

    # Visualize the coefficients...
    for i in range(5,10):
        rank1, rank2 = 4, i
        offset = "0.0,1.0"
        plot_coeff(c, rank1, rank2, labels = labels, offset = offset, rot = 90)
