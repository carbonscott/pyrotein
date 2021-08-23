#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def constant_aminoacid_code():
    aa_dict = {
        "R" : "ARG", "H" : "HIS", "K" : "LYS", "D" : "ASP", "E" : "GLU",
        "S" : "SER", "T" : "THR", "N" : "ASN", "Q" : "GLN", "C" : "CYS",
        "G" : "GLY", "P" : "PRO", "A" : "ALA", "V" : "VAL", "I" : "ILE",
        "L" : "LEU", "M" : "MET", "F" : "PHE", "Y" : "TYR", "W" : "TRP",

        "-" : "MAR"
    }

    return aa_dict




def constant_atomlabel():
    # MAR stands for missing-a-residue;
    # We consider MAR still has 4 placeholder atoms that form a backbone
    label_dict = {
        "ARG" : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
        "HIS" : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
        "LYS" : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'],
        "ASP" : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2'],
        "GLU" : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2'],
        "SER" : ['N', 'CA', 'C', 'O', 'CB', 'OG'],
        "THR" : ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2'],
        "ASN" : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2'],
        "GLN" : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2'],
        "CYS" : ['N', 'CA', 'C', 'O', 'CB', 'SG'],
        "GLY" : ['N', 'CA', 'C', 'O'],
        "PRO" : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD'],
        "ALA" : ['N', 'CA', 'C', 'O', 'CB'],
        "VAL" : ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2'],
        "ILE" : ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1'],
        "LEU" : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2'],
        "MET" : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE'],
        "PHE" : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        "TYR" : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
        "TRP" : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],

        "MAR" : ['N', 'CA', 'C', 'O'],
    }

    return label_dict
