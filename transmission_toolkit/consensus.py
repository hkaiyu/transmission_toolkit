"""Module containing functions for creating consensus sequences"""

# Standard library import
import os
import glob

# Local import
from .parsers import extract_lfv

def get_ref(reference):
    """
    Returns string representation of reference genome given a reference FASTA file.
    If FASTA file has more than one genome in it, returns first genome.
    """    
    consensus = str()

    with open(reference, 'r') as f:
        stop = False
        for line in f.readlines():
            line = line.split('\n')[0]
            if stop and line[0] == '>':
                break
            elif line[0] == '>':
                stop = True
            else:
                consensus += line

    return consensus

def build_consensus(vcf_file, reference):
    """
    With a VCF file and reference file specified, builds consensus sequence.
    """
    consensus = get_ref(reference)
    variants = extract_lfv(vcf_file)
    highest_variants = {}

    for pos in variants:
        highest_freq = 0
        kept_var = None
        for var in pos:
            if variants[pos][var][0] > highest_freq:
                highest_freq = variants[pos][var][0]
                kept_var = var
        highest_variants[pos] = str(kept_var)

    for pos in highest_variants:
        consensus[pos - 1] = highest_variants[pos]
    
    return consensus

def map_consensus(vcfpath, reference):
    """
    Returns a dictionary mapping filenames to its respective consensus.

    Args:
        vcfpath (str): Path to folder containing vcf files
        ref_seq (str): Path to reference genome file
    """
    consensus_dict = dict()
    for vcf in glob.glob(os.path.join(vcfpath,"*.vcf")):
        consensus_dict[vcf] = build_consensus(os.path.join(vcfpath, vcf), reference)
    return consensus_dict
