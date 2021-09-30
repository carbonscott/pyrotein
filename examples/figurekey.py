#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import colorsimple as cs

entry_dict = {
    "A (Rhodopsin)"           : { "shape" : 7, "ps" : 1.5,   "clr" : "#FFB8B8" },
    "B1 (Secretin)"           : { "shape" : 7, "ps" : 1.5,   "clr" : "#00A600" },
    "C (Glutamate)"           : { "shape" : 7, "ps" : 1.5,   "clr" : "#0080FF" },
    "F (Frizzled)"            : { "shape" : 7, "ps" : 1.5,   "clr" : "#AB00FF" },
    "Inactive"                : { "shape" : 6, "ps" : 1.5,   "clr" : "black" },
    "Intermediate"            : { "shape" : 8, "ps" : 1.5,   "clr" : "black" },
    "Active"                  : { "shape" : 2, "ps" : 1.5,   "clr" : "black" },
    "Resolution (>3.5 {\305})": { "shape" : 7, "ps" : 0.5, "clr" : "black" },
}


color_dict = {}
for k, v in entry_dict.items():
    title = k
    shape = v["shape"]
    ps    = v["ps"]
    clr   = v["clr"]

    color_dict[title] = {
        "style" : f"u 1:2 w p pt {shape} ps {ps} lw 1.0 lc rgb '{clr}' title '{title}'",
        "entry" : [],
    }


cs.color_table(color_dict, filename = "xfam-loop.color_table")
