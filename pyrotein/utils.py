#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


def bin_image(img_orig, bin = 4, mode = 1, nan_replace = 0):
    ''' Bin an image for faster display.
    '''
    Y, X = img_orig.shape

    if mode == 0:
        img_bin = []
        for i in range(0, Y, bin):
            for j in range(0, X, bin):
                sub_img = img_orig[i : min(i + bin, Y), j : min(j + bin, X)]
                if np.all(np.isnan(sub_img)):
                    img_bin(append( (i, j, 0) ))
                else:
                    img_bin.append( (i, j, np.nanmean(sub_img)) )

    if mode == 1:
        img_bin = []
        for i in range(0, Y, bin):
            img_bin_y = []
            for j in range(0, X, bin):
                sub_img = img_orig[i : min(i + bin, Y), j : min(j + bin, X)]
                if np.all(np.isnan(sub_img)): 
                    img_bin_y.append( nan_replace )
                else:
                    img_bin_y.append( np.nanmean(sub_img) )
            img_bin.append(img_bin_y)

    return np.array(img_bin)


def read_file(file, numerical = False):
    '''Return all lines in the user supplied parameter file without comments.
    '''
    lines = []
    with open(file,'r') as fh:
        for line in fh.readlines():
            # Separate entries by spaces and remove commented lines...
            words = line.replace('#', ' # ').split()

            # Omit any thing coming after the pound sign in a line...
            if "#" in words: words = words[  : words.index("#")]

            # Save non-empty line...
            if numerical: words = [ float(word) for word in words ]
            if len(words) > 0: lines.append(words)

    return lines
