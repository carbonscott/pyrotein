#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from .constants  import constant_atomlabel, constant_aminoacid_code
from .utils import chunker, sort_dict_by_key
from difflib import ndiff

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




def readonly(file):
    ''' Extract atomic information for every atom and save them in tuples.  
        Users are not encouraged to change values in the results returned.  
    '''
    lines = []
    with open(file,'r') as fh:
        for line in fh.readlines():
            # Skip lines not starting with ATOM or HETATM...
            rec_type = line[0:6].strip()
            if not rec_type in "ATOM HETATM".split(): continue

            # Split a line based on columns defined according to PDB format...
            lines.append( tuple(split(line)) )

    return tuple(lines)




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
        chain_dict =  create_lookup_table(atoms_pdb)
        chain_dict["A"][1002]["CA"]
        ```
    '''
    chain_dict = {}
    for a in atoms_pdb:
        # Decouple parameters...
        resi  = a[6]
        name  = a[2]
        chain = a[5]

        # Initialize by chain...
        if not chain in chain_dict: chain_dict[chain] = {}

        # Initialize by resi for a chain... 
        if not resi  in chain_dict[chain]: chain_dict[chain][resi] = {}

        # Store by exclusive atom name...
        # Atom names are not allowed to have duplicate items
        chain_dict[chain][resi][name] = a

    return chain_dict




# [[[ Dictionary based methods ]]]

def filter_by_resi(chain_dict, nterm, cterm):
    ''' Filter out a subset of an amino acid chain in a range specified by two
        numbers.  

        Use chain_dict = atom_dict[chain] prior to use this function.  
    '''
    return { k : v for k, v in chain_dict.items() if nterm <= k <= cterm }




def filter_by_resn(chain_dict, resn):
    ''' Filter out a subset of an amino acid chain by the residue name
        specified in `resn`.  

        Use chain_dict = atom_dict[chain] prior to use this function.  
    '''
    return { k : v for k, v in chain_dict.items() if "CA" in v and v["CA"][4] == resn }




def resn_to_resi(chain_dict):
    ''' Extract resn to resi mapping (like the sequence viewer on PyMol)
    '''
    aa_dict = { v : k for k, v in constant_aminoacid_code().items() }

    resn_to_resi_dict = {}
    for resi, atom_dict in chain_dict.items():
        atom_first = list(atom_dict.keys())[0]
        resn = atom_dict[atom_first][4]

        if not resn in aa_dict: continue
        resn_to_resi_dict[resi] = aa_dict[resn]

    resn_to_resi_dict = sort_dict_by_key(resn_to_resi_dict)
    return resn_to_resi_dict




def seqi_to_resi(chain_dict, tar_seq, nseqi, cseqi):
    ''' Map seqi to resi for the sequence tar_seq bound by nseqi and cseqi.

        resi is infered by resn_to_resi_dict.  

        Key step is to recognize the lower bound resi that corresponds to
        nseqi.
    '''
    # Initialize mapping...
    seqi_to_resi_dict = { k : None for k in range(nseqi, cseqi + 1) }

    # Extract resn to resi mapping (like the sequence viewer on PyMol)...
    # Non amino acid residue (like ligand) are bypassed
    resn_to_resi_dict = resn_to_resi(chain_dict)

    # Select the bound sequence by nseqi and cseqi...
    tar_seq_bound           = tar_seq[nseqi : cseqi + 1]
    tar_seq_bound_continous = tar_seq[nseqi : cseqi + 1].replace('-', '')

    # Obtain the original sequence from PDB...
    # May have overhead in the beginning or the end
    seq_orig = ''.join([ v for v in resn_to_resi_dict.values() ])

    # Obtain the starting index by string match...
    lb_term = seq_orig.find(tar_seq_bound_continous)

    # Warn user about illegal input sequence...
    if lb_term == -1: 
        print(f"Illegal input sequence!!!")
        print(f"-------------------------")
        for s1, s2 in zip(chunker(tar_seq_bound_continous, 60), chunker(seq_orig, 60)):
            for di, dv in enumerate(ndiff([s1], [s2])): print(dv)
        return seqi_to_resi_dict

    # Obtain the ending index by the length of the coutinous (no '-') sequence...
    ub_term = lb_term + len(tar_seq_bound_continous)

    # Obtain list of resi bound by nseqi and cseqi...
    resi_list       = [ v for v in resn_to_resi_dict.keys() ]
    resi_bound_list = resi_list[lb_term : ub_term]

    # Counter to go through the bound sequence by nseqi and cseqi...
    res_counter = 0

    # Loop through
    for i, seqi in enumerate(range(nseqi, cseqi + 1)):
        # Skip the '-' residue...
        if tar_seq_bound[i] == '-': continue

        # Access the resi...
        resi = resi_bound_list[res_counter]

        # Record the mapping...
        seqi_to_resi_dict[seqi] = resi

        # Increment the residue counter...
        res_counter += 1

    return seqi_to_resi_dict




def extract_xyz_by_seq(atoms_to_extract, chain_dict, tar_seq, nseqi, cseqi):
    ''' Use super_seq as framework to extract coordinates.  Each tar_seq is 
        considered as a subset of super_seq.  

        Sequence alignment directly determines the structure distance matrix.  
    '''
    # Obtain the seq to resi mapping...
    seqi_to_resi_dict = seqi_to_resi(chain_dict, tar_seq, nseqi, cseqi)

    # Obtain size of the seqstr...
    len_chain = cseqi - nseqi + 1

    # Obtain chain info...
    len_res = len(atoms_to_extract)
    len_seq = len_chain * len_res

    # Preallocate memory for storing coordinates...
    xyzs    = np.zeros((len_seq, 3))    # Initialize coordinate matrix
    xyzs[:] = np.nan                    # np.nan for any missing residue

    # Go through each seqi...
    for i, (seqi, resi) in enumerate(seqi_to_resi_dict.items()):
        # From each atom
        for j, atm in enumerate(atoms_to_extract):
            # Derive the matrix index...
            mat_i = i * len_res + j

            # Assign coordinates to matrix at index mat_i...
            if resi in chain_dict:
                if atm in chain_dict[resi]:
                    # 8:8+3 => x,y,z
                    xyzs[mat_i] = chain_dict[resi][atm][8:8+3]

    return xyzs




def extract_xyz_by_atom(atoms_to_extract, chain_dict, nterm, cterm):
    ''' Extract atomic coordinates of interest (specified in the first
        argument) in the lookup table format.  

        if atoms_to_extract is empty `[]`, then it is derived from the resn.  

        Use chain_dict = atom_dict[chain] prior to use this function.  

        All entry must have the same length for RMSD and SVD analysis.
    '''
    # Define atoms used for distance matrix analysis...
    len_backbone = (cterm - nterm + 1) * len(atoms_to_extract)

    # Preallocate memory for storing coordinates...
    xyzs    = np.zeros((len_backbone, 3))    # Initialize coordinate matrix
    xyzs[:] = np.nan                         # np.nan for any missing residue

    # From each residue
    for i in range(nterm, cterm + 1):
        # From each atom
        for j, p in enumerate(atoms_to_extract):
            # Derive the matrix index...
            mat_i = (i - nterm) * len(atoms_to_extract) + j

            # Assign coordinates to matrix at index mat_i...
            if i in chain_dict:
                if p in chain_dict[i]:
                    # 8:8+3 => x,y,z
                    xyzs[mat_i] = chain_dict[i][p][8:8+3]

    return xyzs




def extract_xyz_with_sidechain(tar_seq, super_seq, chain_dict, nterm, cterm):
    ''' Extract xyzs from a protein chain, whose sequecne is tar_seq.  Two
        sequecnes are considered in this function.  super_seq is used to 
        distinguish three scenarios about how to extract xyzs.  

        - missing : - != Q (tar != super)
        - mismatch: N != Q (tar != super)
        - match   : Q == Q (tar == super)

        In the mismatch scenario, the resn in tar_seq only contributes to main
        chain atoms (N, CA, C, O) if they exit.  

        In the missing scenario, no xyz is extracted from tar_seq.

        Use chain_dict = atom_dict[chain] prior to use this function.  
    '''
    # Import label_dict and aa_dict...
    label_dict = constant_atomlabel()
    aa_dict    = constant_aminoacid_code()

    # Count atoms used for distance matrix analysis based on super_seq...
    len_xyzs = np.sum( [ len(label_dict[aa_dict[i]]) for i in super_seq ] )

    # Preallocate memory for storing coordinates...
    xyzs    = np.zeros((len_xyzs, 3))    # Initialize coordinate matrix
    xyzs[:] = np.nan                     # np.nan for any missing residue

    # From each residue
    mat_i = 0
    for seqi, resi in enumerate(range(nterm, cterm + 1)):
        # Extract the resn in both sequences...
        tar_atoms   = label_dict[aa_dict[  tar_seq[seqi]]]
        super_atoms = label_dict[aa_dict[super_seq[seqi]]]

        # Distingusih scenarios
        # If missing, skip the residue???
        if tar_seq[seqi] == '-':
            mat_i += len(super_atoms)
            continue

        # If mismatch, consider main chain and skip???
        len_skip = 0
        if tar_seq[seqi] != super_seq[seqi]:
            # Consider mainchain atoms only...
            tar_atoms = label_dict["GLY"]

            # Get the skip length...
            len_skip = len(super_atoms) - len(tar_atoms)

            ## print(f"{seqi}|{resi} : {tar_seq[seqi]}|{super_seq[seqi]} : {len(tar_atoms)}|{len(super_atoms)}")

        # If match???
        resi_dict = chain_dict[resi]
        for j, atm in enumerate(tar_atoms):
            # 8:8+3 => x,y,z
            if atm in resi_dict:
                xyzs[mat_i] = resi_dict[atm][8:8+3]
            mat_i += 1
        mat_i += len_skip

    return xyzs




def is_standard(atom1_head, atom1_tail, atom2_head, atom2_tail):
    ''' Check if two pose vectors vector1 and vector2 form an acute angle (< 90 degree).  
        If not, swap head and tail in vector1.  

        4 args should be numpy arrays.  
    '''
    # Set the default result to be True...
    res = True

    # Construct vectors...
    vec1 = atom1_head - atom1_tail
    vec2 = atom2_head - atom2_tail

    # False if not standard...
    if np.dot(vec1, vec2) < 0: res = False

    return res




def standardize_sidechain(atom_dict):
    ''' Standardize side chain atoms for a given atomic structure in atom_dict.
        It facilitates distance matrix analysis of sidechains.  The order of atoms 
        should follow a standard.  The type of atoms don't make a difference in 
        distance matrix analysis.  
    '''
    # Consider resiudes with ambiguous atomic position...
    # e.g. NH1 and NH2 in ARG can be swapped
    ambi_dict = {
        "ARG" : [ "CD", "NE", "NH1", "NH2" ],
        "ASP" : [ "CA", "CB", "OD1", "OD2" ],
        "ASN" : [ "CA", "CB", "OD1", "ND2" ],
        "GLU" : [ "CB", "CG", "OE1", "OE2" ],
        "GLN" : [ "CB", "CG", "OE1", "NE2" ],
        ## "HIS" : [ "CA", "CB", "ND1", "CD2" ],
        ## "THR" : [ "C" , "CA", "OG1", "CG2" ],
        "VAL" : [ "C" , "CA", "CG1", "CG2" ],
        "LEU" : [ "CA", "CB", "CD1", "CD2" ],
        "PHE" : [ "CA", "CB", "CD1", "CD2" ],
        "TYR" : [ "CA", "CB", "CD1", "CD2" ],
    }

    # Specify the swapping rule to standardize sidechains...
    swap_dict = {
        "ARG" : [["NH1", "NH2"]],
        "ASP" : [["OD1", "OD2"]],
        "ASN" : [["OD1", "ND2"]],
        "GLU" : [["OE1", "OE2"]],
        "GLN" : [["OE1", "NE2"]],
        ## "HIS" : [["ND1", "CD2"], ["CE1", "NE2"]],
        ## "THR" : [["OG1", "CG2"]],
        "VAL" : [["CG1", "CG2"]],
        "LEU" : [["CD1", "CD2"]],
        "PHE" : [["CD1", "CD2"], ["CE1", "CE2"]],
        "TYR" : [["CD1", "CD2"], ["CE1", "CE2"]],
    }

    # Fix the ordering throughout atom_dict...
    for chain, chain_dict in atom_dict.items():
        for resi, resi_dict in chain_dict.items():
            # Skip entries not containing CA backbone...
            if not "CA" in resi_dict: continue

            # Fetch the name of residue...
            # 4 => resname, check `pr.atom.spec()`
            resn = resi_dict["CA"][4]

            # Skip entries not having ambiguous atom placement...
            if not resn in ambi_dict: continue

            # Skip entries that have missing atoms...
            if not np.all([ i in resi_dict for i in ambi_dict[resn] ]): continue

            # Get pose vectors...
            # 8:8+3 => x, y, z
            atom1_head, atom1_tail, atom2_head, atom2_tail = \
            [ np.array(resi_dict[i][8:8+3]) for i in ambi_dict[resn] ]

            # Swap xyz when pose is not standard...
            if not is_standard(atom1_head, atom1_tail, atom2_head, atom2_tail):
                for atom1_swap, atom2_swap in swap_dict[resn]:
                    resi_dict[atom1_swap][8:8+3], resi_dict[atom2_swap][8:8+3] = \
                    resi_dict[atom2_swap][8:8+3], resi_dict[atom1_swap][8:8+3]

    return None




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
    resi_dict = extract_residue(atoms)

    # Fill the table (atoms_track_table)...
    for resi in resi_dict.keys():
        for atom_to_track in atoms_to_track:
            if resi not in atoms_track_table[atom_to_track]: 
                atoms_track_table[atom_to_track][resi] = 0

    return atoms_track_table




def get_coord(atom_select): return atom_select[8:8+3]




