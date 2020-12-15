#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

''' The program extracts atomic information from a truncated pdb file that is
    exported from PyMol after a selection statement is executed, and thus has a
    format below (part of the file is shown).

HETATM    1  CAA CYC A 175      31.603  38.077  18.635  1.00 18.22           C 
HETATM    2  CAB CYC A 175      23.511  33.969  15.096  1.00 13.95           C 
HETATM    3  CAC CYC A 175      33.943  30.709  11.039  1.00 22.02           C 
HETATM    4  CAD CYC A 175      36.012  36.289  16.517  1.00 23.57           C 
HETATM    5  NA  CYC A 175      30.590  35.666  15.935  1.00 17.82           N 
HETATM    6  CBA CYC A 175      32.240  37.352  19.790  1.00 21.28           C 
HETATM    7  CBB CYC A 175      22.926  32.891  15.989  1.00 13.87           C 
HETATM    8  CBC CYC A 175      33.396  29.891   9.878  1.00 17.93           C 
HETATM    9  CBD CYC A 175      36.524  37.532  15.828  1.00 27.36           C 
HETATM   10  NB  CYC A 175      26.189  35.771  16.906  1.00 13.89           N
'''

def split(line):
    ''' It breaks a sing string into a list of pdb coordinates.  

        Protein data bank format:
        https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    '''
    return [ line[0 : 6].strip(),    #  0  ATOM or HETATM
             line[6 :12].strip(),    #  1  Atom serial number (1..<whatever>)
             line[12:16].strip(),    #  2  Atom name (CA)
             line[16   ].strip(),    #  3  Alternate location indicator
             line[17:21].strip(),    #  4  Residue name
             line[21   ].strip(),    #  5  Chain identifier

             int  (line[22:26]),     #  6  Residue sequence number
             line[26   ].strip(),    #  7  Code for insertions of residues
             float(line[30:38]),     #  8  X
             float(line[38:46]),     #  9  Y
             float(line[46:54]),     # 10  Z
             float(line[54:60]),     # 11  Occupancy
             float(line[60:66]),     # 12  Temperature factor

             line[72:76].strip(),    # 13  Segment identifier
             line[76:78].strip(),    # 14  Element symbol 
           ]



def read(file):
    ''' Extract atomic information for every atom.  
    '''
    lines = []   

    with open(file,'r') as fh:
        for line in fh.readlines():
            # Skip lines not starting with ATOM or HETATM...
            rec_type = line[0:6].strip()
            if not rec_type in "ATOM HETATM".split(): continue

            # Split a line based on columns defined according to PDB format...
            lines.append( split(line) )

    return lines




def spec():
    ''' Print out the PDB format formating rule.
    '''
    print( '''
    The PDB formating rule: 
    _______________________ 

           0  ATOM or HETATM
           1  Atom serial number (1..<whatever>)
           2  Atom name (CA)
           3  Alternate location indicator
           4  Residue name
           5  Chain identifier

           6  Residue sequence number
           7  Code for insertions of residues
           8  X
           9  Y
          10  Z
          11  Occupancy
          12  Temperature factor

          13  Segment identifier
          14  Element symbol 
           '''
         )




def create_lookup_table(atoms_pdb):
    ''' Create a lookup table to access atomic coordinates.  

        ```
        atom_dict =  create_lookup_table(atoms_pdb)
        atom_dict["A"][1002]["CA"]
        ```
    '''
    atom_dict = {}
    for a in atoms_pdb:
        # Decouple parameters...
        resi  = a[6]
        name  = a[2]
        chain = a[5]

        # Initialize by chain...
        if not chain in atom_dict: atom_dict[chain] = {}

        # Initialize by resi for a chain... 
        if not resi  in atom_dict[chain]: atom_dict[chain][resi] = {}

        # Store by exclusive atom name...
        # Atom names are not allowed to have duplicate items
        atom_dict[chain][resi][name] = a

    return atom_dict




# [[[ Dictionary based methods ]]]

def extract_segment(atom_dict, chain, nterm, cterm):
    ''' Extract a subset of an amino acid chain in a range specified by two
        numbers.  

        Returns atmoic coordinates in a lookup table format.  
    '''
    return { k : v for k, v in atom_dict[chain].items() if nterm <= k <= cterm }




def extract_backbone_xyz(atom_dict, chain, nterm, cterm):
    ''' Extract the backbone atoms ["N", "CA", "C", "O"] from atomic
        coordinates in the lookup table format.  
    '''
    # Extract the segment of amino acids...
    rho_dict = extract_segment(atom_dict, chain, nterm, cterm)

    # Define atoms used for distance matrix analysis...
    backbone = ["N", "CA", "C", "O"]
    len_backbone = (cterm - nterm + 1) * len(backbone)

    # Preallocate memory for storing coordinates...
    xyzs    = np.zeros((len_backbone, 3))    # Initialize coordinate matrix
    xyzs[:] = np.nan                         # np.nan for any missing residue

    # From each residue
    for i in range(nterm, cterm + 1):
        # From each backbone atom
        for j, p in enumerate(backbone):
            # Derive the matrix index...
            mat_i = (i - 1) * len(backbone) + j

            # Assign coordinates to matrix at index mat_i...
            if i in rho_dict:
                if p in rho_dict[i]:
                    xyzs[mat_i] = get_coord(rho_dict[i][p])

    return xyzs




# [[[ List based methods ]]]

def import_connect(fl_pdb):
    ''' Extract connection information for a molecule.  
    '''
    lines = []   

    with open(fl_pdb,'r') as fh:
        for line in fh.readlines():
            # Skip lines not starting with ATOM or HETATM...
            rec_type = line[0:7].strip()
            if not rec_type in "CONECT".split(): continue

            # Create an space-filled string with 32 characters long...
            s = line + " " * ( 32 - len(line) )

            # Split a line based on columns defined according to PDB format...
            lines.append( [
                line[0:7].strip(),
                line[7:12].strip(),
                line[12:17].strip(),
                line[17:22].strip(),
                line[22:27].strip(),
                line[27:32].strip(),
            ] )

    return lines




def get_chain(atoms, chain):
    ''' Return a list of atoms from a user-defined chain
    '''
    return [ atom for atom in atoms if atom[0] == "ATOM" and atom[5] == chain ]




def extract_residue(atoms):
    ''' Convert a list of atoms to a look-up table for amino acid residues.

        Return { 20 : "LYS", 21 : "PRO", ... }.
    '''
    resi_dict = {}
    for atom in atoms:
        resi_name = atom[4]
        resi_num  = atom[6]
        if resi_num not in resi_dict: resi_dict[resi_num] =  resi_name

    return resi_dict




def track_atoms(atoms, atoms_to_track):
    ''' Return a table of atoms to track for each residue.

            0 1 2 3 <-- residue numbers
        CA  1 1 1 1
        C   0 1 1 1
        O   1 1 1 0
        N   1 1 1 1
        ^
        |__ atoms to track

        - 1 stands for existence;
        - 0 stands for absence;
    '''
    # Initialize tracking table...
    atoms_track_table = {}
    for atom_to_track in atoms_to_track:
        atoms_track_table[atom_to_track] = \
            { atom[6] : 1 for atom in atoms if atom[2] == atom_to_track }

    # Obtain all residues...
    resi_dict = get_resname(atoms)

    # Fill the table (atoms_track_table)...
    for resi in resi_dict.keys():
        for atom_to_track in atoms_to_track:
            if resi not in atoms_track_table[atom_to_track]: 
                atoms_track_table[atom_to_track][resi] = 0

    return atoms_track_table




def get_coord(atom_select): return atom_select[8:8+3]




