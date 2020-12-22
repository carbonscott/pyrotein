#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import openpyxl
import colorsimple as cs


def load_xlsx(fl_input):
    # Load the spreadsheet...
    bk = openpyxl.load_workbook(fl_input, data_only = True)
    st = bk["Sheet1"]

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

        entries.append(val)

    return entries




def label_TMs():
    return {"TM1"  : [ 33,  65],
            "TM2"  : [ 70, 100],
            "TM3"  : [105, 140],
            "TM4"  : [149, 173],
            "TM5"  : [199, 236],
            "TM6"  : [240, 277],
            "TM7"  : [288, 307],
            "H8"   : [310, 322]}




def color_species(items, s = 50, v = 100):
    ''' 
    '''
    assert len(set(items)) == len(items), \
        "Duplicate item is not allowed in the input."

    # Get number of colors...
    num = len(items)

    # Divide the color palette...
    div = int(360 / num)

    # Assign color to each item...
    color_dict = {}
    for i, item in enumerate(items): 
        color_dict[item] = "#" + cs.hsv_to_hex(i * div, s, v)

    return color_dict




def color_table(color_dict):
    import GnuplotPy3
    gp = GnuplotPy3.GnuplotPy3()

    gp( "set terminal postscript eps  size 3.5, 2.62 \\")
    gp( "                             enhanced color \\")
    gp( "                             font 'Helvetica,14' \\")
    gp( "                             linewidth 2")
    gp(f"set output 'color_table.eps'")
    gp("set xrange [1:2]")
    gp("set yrange [1:2]")
    gp("unset border")
    gp("unset xtics")
    gp("unset ytics")

    gp("plot \\")
    for n, c in color_dict.items():
        gp(f"'-' using 1:2 with points pointtype 7 linecolor rgb '{c}' title '{n}',\\")
    gp("")

    for i in range(len(color_dict)):
        gp(f"0 0")
        gp( "e")
    gp("exit")

    return None


## color_dict = color_species(range(7))
## color_table(color_dict)
