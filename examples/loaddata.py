#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import openpyxl


def load_xlsx(fl_input, sheet = "Sheet1"):
    # Load the spreadsheet...
    bk = openpyxl.load_workbook(fl_input, data_only = True)
    st = bk[sheet]

    # Fetch entries from the 2nd row...
    entries = []
    for i in st.iter_rows(min_row = 2):
        # Unpack all values
        val = [ j.value for j in i ]

        # Stop iteration if row starts with None
        if val[0] is None: break

        # Unpack parameters
        select, pdb, chain, species = val[:4]

        # Skip unselected...
        if select == "no": 
            print(f"{pdb}-{chain}-{species} is not selected.")
            continue

        # Go through each chain and record...
        for eachchain in chain.split(): 
            # Replace the chain column with a single chain identifier...
            val[2] = eachchain
            entries.append(val.copy())

    return entries




def label_TMs():
    return {"TM1"  : [ 35,  50],
            "TM2"  : [ 70, 100],
            "TM3"  : [105, 140],
            "TM4"  : [149, 173],
            "TM5"  : [199, 236],
            "TM6"  : [240, 277],
            "TM7"  : [288, 307],
            "H8"   : [310, 322]}

