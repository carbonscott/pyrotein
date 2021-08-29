#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' This program process the following request.  

    Each request takes two input sequence: 
    - one sequence as reference that is typically from structure-based
      sequence alignment.  
        - db_fasta: A file that specifies all reference fasta sequences.  
                    The framework for data correspondance is presumbly 
                    taken from here.  

    - one sequence to align with the reference sequence.  
        - proc_fasta: A file that includes all sequences that are going to 
                      align with its corresponding reference sequence.  

    The reference sequence contains '-' as it is often an output from a
    sequence alignment program.  Convert all '-' into 'X' to trick the service
    called needle, providing pairwise sequence alignment, from
    https://www.ebi.ac.uk/Tools/psa/emboss_needle/ .  

    How to know which the reference sequence is, given any sequence to align?  
    - UniProt name and species name work together as a key to find the 
      reference sequence.  
    - e.g. "opsd_bovin" 
            ~~~~ ^^^^^
              |____|____ UniProt name
                   |____ Species label
    - Special label is a 5-letter code, retrivable according to 'organism.dat' 
      by the common name (C name).  e.g. bovine (C) --> bovin (label)

    The output sequence will be saved in two directories
    - psa.i: temp files as inputs for the needle program
    - psa.o: output files from needle, containing xxx.aln.txt file.  

    This program goes through items marked as "yes" in fl_xlsx, a database by 
    which user can choose if an entry should be selected.  


    Auxiliary inputs:
    - "organism.dat": the UniProt naming convention is fetched here.  
                      Downloaded from https://www.uniprot.org/docs/speclist
                      on Aug 27, 2021.  
    - fl_xlsx ("gpcrdb.all.xlsx"): the database file.  


    Caveats:

    - Use aln.s2.checkdash.py to check if there is '-' in the output reference
      sequence in the xxx.aln.txt file.  If so, the sequence to align has an
      unconventional amino acid members or order.  

    - Prepare a new proc_fasta if new/updated sequecnes are accounted for.  

    - Run aln.s5.remove_defect_fasta.py to remove wrong sequence files if an
      update is demanded.  Otherwise, the program respects existing files and 
      will bypass the process by default.  

    - (Optional) Fetching refined structures from gpcrdb
      - aln.s3.dl_refine.py
      - aln.s4.fetch_fasta.py

'''

import os
import pyrotein as pr
from loaddata import load_gpcrdb_xlsx, load_uniprot_species, name_to_code
import multiprocessing as mp

# Give a name to the analysis...
job_name = "xfam"

# Read uniprot specification...
fl_spec   = "organism.dat"
spec_dict = load_uniprot_species(fl_spec)
name_to_code_dict = name_to_code(spec_dict, "c")    # "c" stands for common name

# Read the db fasta...
db_fasta = f"gpcrdb.{job_name}.fasta"
frame_dict = pr.fasta.read(db_fasta)

# Read the fasta to process...
proc_fasta = "step2.interest.fasta"
proc_dict = pr.fasta.read(proc_fasta)

# Test the concept of fetching sequence alignment 

# Specify chains to process...
fl_xlsx    = f"gpcrdb.all.xlsx"
sheet      = f"{job_name}"
lines      = load_gpcrdb_xlsx(fl_xlsx, sheet = sheet)
drc        = f"pdb"
drc_input  = f"{job_name}.psa.i"
drc_export = f"{job_name}.psa.o"

if False:
    for line in lines[:]:
        # Unpack parameters
        uniprot = line[1].lower()
        spec    = line[5].lower()
        pdb     = line[7].lower()
        chain   = line[10].upper()

        # Form entry name...
        entry   = f"{pdb}_{chain}"

        # Derive the id like uniprot_speccode...
        speccode = name_to_code_dict[f"{spec}"]
        entry_ref = f"{uniprot}_{speccode}"
        print(f"Processing {entry} from {entry_ref}...")

        # Find the reference according to species (uniprot_speccode)...
        fl_ref_fasta = os.path.join(drc_input, f"tmp.ref.{entry_ref}.fasta")

        # Only Write file when needed
        # "-" to "X" is the trick that provides placeholder for residue
        if not entry_ref in frame_dict: 
            print(f"{entry_ref} N/A!!!  Skip this species...")
            continue
        if not os.path.exists(fl_ref_fasta): pr.fasta.write(fl_ref_fasta, entry_ref, frame_dict[entry_ref].replace("-", "X"))

        # Export fasta as input sequence...
        fl_i = os.path.join(drc_input, f"tmp.{entry}.fasta")

        # Only Write file when needed
        if os.path.exists(fl_i): continue
        pr.fasta.write(fl_i, entry, proc_dict[entry])

        # Extract sequence...
        fl_export = os.path.join(drc_export, f"out.{entry}")

        ## print(f"./emboss_needle.py --outfile={fl_export} --asequence={fl_ref_fasta} --bsequence={fl_i} --gapopen=10 --gapext=0.5 --email=cwang230@uic.edu --format=fasta")
        os.system(f"./emboss_needle.py --outfile={fl_export} --asequence={fl_ref_fasta} --bsequence={fl_i} --gapopen=10 --gapext=0.5 --email=cwang230@uic.edu --format=fasta")

if True:
    def psa_run(line):
        # Unpack parameters
        uniprot = line[1].lower()
        spec    = line[5].lower()
        pdb     = line[7].lower()
        chain   = line[10].upper()

        # Form entry name...
        entry   = f"{pdb}_{chain}"

        # Derive the id like uniprot_speccode...
        if not spec in name_to_code_dict: 
            print(f"!!! UniProt entry error -- No {spec}")
            return None
        speccode = name_to_code_dict[spec].strip()
        entry_ref = f"{uniprot}_{speccode}"
        print(f"Processing {entry} from {entry_ref}...")

        # Find the reference according to species (uniprot_speccode)...
        fl_ref_fasta = os.path.join(drc_input, f"tmp.ref.{entry_ref}.fasta")

        # Only Write file when needed
        # "-" to "X" is the trick that provides placeholder for residue
        if not entry_ref in frame_dict: 
            print(f"!!! Reference entry error -- No {entry_ref} ")
            return None
        if not os.path.exists(fl_ref_fasta): pr.fasta.write(fl_ref_fasta, entry_ref, frame_dict[entry_ref].replace("-", "X"))

        # Export fasta as input sequence...
        fl_i = os.path.join(drc_input, f"tmp.{entry}.fasta")

        # Only Write file when needed
        if not entry in proc_dict: 
            print(f"!!! PDB entry error -- No {entry}")
            return None
        if not os.path.exists(fl_i): pr.fasta.write(fl_i, entry, proc_dict[entry])

        # Extract sequence...
        fl_export = os.path.join(drc_export, f"out.{entry}")

        ## print(f"./emboss_needle.py --outfile={fl_export} --asequence={fl_ref_fasta} --bsequence={fl_i} --gapopen=10 --gapext=0.5 --email=cwang230@uic.edu --format=fasta")
        if not os.path.exists(f"{fl_export}.aln.txt"): os.system(f"./emboss_needle.py --outfile={fl_export} --asequence={fl_ref_fasta} --bsequence={fl_i} --gapopen=10 --gapext=0.5 --email=cwang230@uic.edu --format=fasta")

        return None


    if __name__ == "__main__":
        num_job = 4
        with mp.Pool(num_job) as proc:
            proc.map( psa_run, lines[:] )
