#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pymolPy3
import os
from pmview import view_dict

# Specify the file to visualize...
drc      = "pdb.obs"
pdb      = "6wiv"
path_pdb = os.path.join(drc, f"{pdb}.pdb")

# Launch pymol...
pm = pymolPy3.pymolPy3()

# Load presets...
pm("bg white")
pm("set cartoon_fancy_helices, 1")
pm("set cartoon_highlight_color, grey90")
pm("set cartoon_dumbbell_length, 1")
pm("set cartoon_dumbbell_width, 0.3")
pm("set cartoon_dumbbell_radius, 0.2")
pm("set sphere_scale, 0.3")
pm("window size, 1500, 1500")

# Load molecule...
pm(f"load {path_pdb}")

# Set view...
pm(f"{view_dict[pdb]['view']}")

# Set color of protein...
pm(f"set cartoon_color, {view_dict[pdb]['protein_color']}")

# Select ligands...
pm(f"select lig, not polymer.protein")
pm("disable %lig")

# Set color of ligands...
pm(f"set stick_color, {view_dict[pdb]['ligand_color']}")

input()
