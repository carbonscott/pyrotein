#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .utils import get_key_by_max_value


def read(fl_fasta):
    ''' Extract sequence from a fasta file.  
    '''

    seq = {}
    with open(fl_fasta,'r') as fh:
        for line in fh.readlines():
            if line.startswith(">"): 
                k = line[1:].rstrip()
                seq[k] = ""
            else: seq[k] += line.rstrip()

    return seq



def mask_pairseq(seq1, seq2, null = '-'):
    ''' If both seq1 and seq2 have non-null value at a position, the value at
        the position is associated with True.  
    '''
    mask = {}
    for i in range(len(seq1)):
        mask[i] = False if '-' in seq1[i] + seq2[i] else True
    return mask




def create_seqmask(seq, null = '-'):
    ''' Map False to '-' and True to non '-' in seq index.
    '''
    mask = {}
    for i in range(len(seq)):
        mask[i] = False if '-' in seq[i] else True
    return mask




def seq_to_resi(seq, resi_non_null, null = '-'):
    ''' Map seq index to resi.  
        seq is a sequence.  
        resi_non_null is the resi to the first non '-' residue.  
    '''
    seqmask = create_seqmask(seq, null = null)

    id_aux = resi_non_null
    seq_to_resi_dict = {}
    for k, v in seqmask.items():
        if v :
            seq_to_resi_dict[k] = id_aux
            id_aux += 1
        else:
            seq_to_resi_dict[k] = None
    return seq_to_resi_dict




def tally_resn_in_seqs(seq_dict):
    ''' Tally the occurence of each residue from a sequence alignment fasta file.
        The input is a seqeuence dictionary.  
    '''
    tally_dict = {}

    for k, v in seq_dict.items():
        for i, resi in enumerate(v):
            # Initialize at resi position at i...
            if not i in tally_dict: tally_dict[i] = {}

            # Count 1 when resi was found the first time...
            if not resi in tally_dict[i]: tally_dict[i][resi] = 1
            else: tally_dict[i][resi] += 1

    return tally_dict




def infer_super_seq(tally_dict):
    ''' Infer the most representative residue based on a tallied result (dict).  
    '''
    return ''.join( [ get_key_by_max_value(v) for v in tally_dict.values() ] )




def read_constant_aminoacid_code():
    aa_dict = {
        "R" : "ARG", "H" : "HIS", "K" : "LYS", "D" : "ASP", "E" : "GLU",
        "S" : "SER", "T" : "THR", "N" : "ASN", "Q" : "GLN", "C" : "CYS",
        "G" : "GLY", "P" : "PRO", "A" : "ALA", "V" : "VAL", "I" : "ILE",
        "L" : "LEU", "M" : "MET", "F" : "PHE", "Y" : "TYR", "W" : "TRP",

        "-" : "MAR"
    }

    return aa_dict




def find_mismatch(ref, tar):
    ''' Return a list of index points to mismatched residue.  Both residues 
        should be available.  That is to say,

        ref and tar must have the same length.
    '''
    len_seq = len(ref)

    mismatch_list = []
    for i in range(len_seq):
        if ref[i] + tar[i] == "--": continue
        if '-' in ref[i] + tar[i]: continue

        if ref[i] != tar[i]: mismatch_list.append(i)

    return mismatch_list
