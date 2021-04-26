#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .utils import get_key_by_max_value
from .atom  import constant_atomlabel, constant_aminoacid_code


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




def seqi_non_null(seq, null = '-'):
    ''' Return a list of sequence index when resn is not '-'.
    '''
    return [ k for k, v in mask_seq(seq, null = null).items() if v ]




def seqi_null(seq, null = '-'):
    ''' Return a list of sequence index when resn is not '-'.
    '''
    return [ k for k, v in mask_seq(seq, null = null).items() if not v ]




def resi_to_seqi(resi, super_seq, nterm):
    ''' Convert resi to seqi.
    '''
    label_dict = constant_atomlabel()
    aa_dict = constant_aminoacid_code()

    return sum( len(label_dict[aa_dict[super_seq[seqi]]]) for seqi in range(resi - nterm) )




def seq_to_resi(seq, resi_non_null, null = '-'):
    ''' Map seq index to resi.  
        seq is a sequence.  
        resi_non_null is the resi to the first non '-' residue.  
    '''
    seqmask = mask_seq(seq, null = null)

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




def diff_seq(tar, ref):
    ''' Return a differnce sequence, in which identical resn will be replaced
        with '-', and different resn will be replaced with the resn in tar.  
    '''
    len_ref = len(ref)
    len_tar = len(tar)

    assert len_ref == len_tar, f"Error in length: len(ref) = {len_ref}; len(tar) = {len_tar}"

    seqdiff = ''
    for i in range(len_ref):
        if ref[i] == tar[i]: seqdiff += '-'
        else:                seqdiff += tar[i]

    return seqdiff




def strip_null(seq, null = '-'): 
    ''' Removing leading and trailing '-' in a sequence.  
    '''
    return seq.strip(null)




def get_lseqi(seq):
    ''' Return the leftmost sequence index corresponds to first non-null ('-') resn.
    '''
    return seq.find( strip_null(seq)[0] )




def get_rseqi(seq):
    ''' Return the rightmost sequence index corresponds to first non-null ('-') resn.
    '''
    return seq.rfind( strip_null(seq)[-1] )
