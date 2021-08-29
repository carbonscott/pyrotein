#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from zipfile import ZipFile
import os
import pyrotein as pr
import requests

# Give a name to the analysis...
job_name = "xfam"

drc = f"{job_name}.pdb.refine"
fl_req_refine = f"{job_name}.req_refine.dat"

# Fetch all PDBs required to refine...
pdbs = [ i[0] for i in pr.utils.read_file(fl_req_refine) ]

# Go through each
for pdb in pdbs[:]:
    # [[[ DOWNLOAD CONTENT ]]]
    urls = [ f"https://gpcrdb.org/structure/refined/{pdb}-complex/download_pdb", 
             f"https://gpcrdb.org/structure/refined/{pdb}_full/download_pdb", ]

    fl_o = os.path.join(drc, f"{pdb}.zip")

    print()
    print(f"Processing {pdb}...")

    for url in urls[:]:
        # Get a response...
        if os.path.exists(fl_o): 
            print(f"~~~ Skipping {url}")
            break
        r = requests.get(url, allow_redirects=True)

        # Try the next available url...
        if not r.status_code == 200: 
            print(f"!!! NA {r.url}")
            r.close()
            continue

        # Write contents...
        open(fl_o, 'wb').write(r.content)

        r.close()

        break



    # [[[ EXTRACT CONTENT ]]]
    fl_pdb = os.path.join(drc, f"{pdb}.pdb")
    if os.path.exists(fl_pdb): 
        print(f"### already there {pdb}")
        continue

    fl_o = os.path.join(drc, f"{pdb}.zip")
    if not os.path.exists(fl_o): 
        print(f"!!! no refined found {pdb}")
        continue

    with ZipFile(f'{fl_o}', "r") as zf:
        fl_list = zf.namelist()
        for fl_each in fl_list:
            if not fl_each.endswith(".pdb"): continue

            # Extract the file to drc...
            zf.extract(fl_each, drc)

            # Rename it right after...
            path_each = os.path.join(drc, fl_each)
            if not os.path.exists(path_each): continue

            newpath_each = os.path.join(drc, f"{pdb.lower()}.pdb")
            os.rename(path_each, newpath_each)

