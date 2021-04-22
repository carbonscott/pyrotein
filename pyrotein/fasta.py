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




def mask_seq(seq, null = '-'):
    ''' Map False to '-' and True to non '-' in seq index.
    '''
    mask = {}
    for i in range(len(seq)):
        mask[i] = False if '-' in seq[i] else True
    return mask




def seq_to_resi(seq_mask, resi_tar):
    ''' Map seq index to resi index.  
    '''
    id_aux = 0
    seq_to_resi_dict = {}
    for k, v in seq_mask.items():
        if v :
            seq_to_resi_dict[k] = resi_tar[id_aux]
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
