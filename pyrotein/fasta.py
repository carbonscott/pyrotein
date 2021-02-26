#!/usr/bin/env python3
# -*- coding: utf-8 -*-


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
