#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pyrotein as pr
from loaddata import load_xlsx, label_TMs
from display import plot_dmat

def get_neighbor(data, baseline, pct_tol):
    ''' Return all residues that satisfies the condition:

        baseline * (1-pct_tol) <= data <= baseline * (1+pct_tol)
    '''
    # Convert the number to percet...
    pct_tol /= 100

    # Find the bounds...
    lower_bound = baseline * (1 - pct_tol)
    upper_bound = baseline * (1 + pct_tol)

    # Obtain residues satisfying conditions...
    items = np.where( (lower_bound <= data) & (data <= upper_bound) )[0]

    return items


def get_lower(data, baseline):
    ''' Return all residues that satisfies the condition:

        data <= baseline
    '''
    # Find the bounds...
    upper_bound = baseline

    # Obtain residues satisfying conditions...
    items = np.where( data <= upper_bound )[0]

    return items


def groupSequence(lst): 
    res = [[lst[0]]] 

    for i in range(1, len(lst)): 
        if lst[i-1]+1 == lst[i]: 
            res[-1].append(lst[i]) 
        else:
            res.append([lst[i]]) 
    return res 


def twoends(groupseq):
    res = []
    for i in groupseq:
        if len(i) > 1: res.append([i[0], i[1]])
    return res


# load the data...
rmsd_dmat = np.load("rmsd_dmat.npy")

# Calculate the mean...
column_mean_dmat = np.nanmean(rmsd_dmat, axis = 0, keepdims = False)
fl_output = "rmsd.column_mean.dat"
with open(fl_output,'w') as fh:
    for i, v in enumerate(column_mean_dmat):
        fh.write(f"{i} {v}")
        fh.write("\n")

# Obtain the min and max...
mean_min, mean_max = np.min(column_mean_dmat), np.max(column_mean_dmat)

# Scan through all mean values and record the residues in the neighborhood...
pct_tol = 1
scan_rate = 1e-2
scan_range = np.arange(mean_min, mean_max + scan_rate, scan_rate)
res_cnt = []
for baseline in scan_range:
    items = get_neighbor(column_mean_dmat, baseline, pct_tol)
    res_cnt.append( [baseline, len(items)] )

fl_output = "rmsd.scan.dat"
with open(fl_output,'w') as fh:
    for x, y in res_cnt:
        fh.write(f"{x} {y}")
        fh.write("\n")
