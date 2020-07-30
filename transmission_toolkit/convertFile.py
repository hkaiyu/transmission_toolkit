"""
Module for converting file types using BioPython's SeqIO module
"""
import os
from Bio import AlignIO
from transmission_toolkit.utils import getpathleaf
from transmission_toolkit.VCFtools import build_majority_consensus, build_minor_consensus

LINE_WRAP = 80 # Max. length of each line in fasta file 

def vcf2fasta(
    vcf_path, 
    reference, 
    output_dir="", 
    line_length=LINE_WRAP, 
    min_AF=0,
    max_AF=1,
    parse_type='biallelic',
    masks=None,
    mask_status='hide',
    consensus_type='majority'
    ):
    """
    Writes a FASTA file with a VCF file and reference FASTA file as input.
    """
    # If user specifies output directory
    if output_dir:

        # Check if directory already exists
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

    name = getpathleaf(vcf_path).split('.')[0]
    path = os.path.join(output_dir, name + '.fna')
    if consensus_type == 'majority':
        seq = build_majority_consensus(
            vcf_path, 
            reference, 
            masks=None, 
            mask_status='hide'
        )
    elif consensus_type == 'minor':
        seq = build_minor_consensus(
            vcf_path, 
            reference, 
            min_AF=0, 
            max_AF=1, 
            masks=None, 
            mask_status='hide'
        )
    else:
        raise ValueError(f'Unexpected consensus_type: {consensus_type}.')

    with open(path, 'w') as f:
        f.write(f'>{name}\n')
        for i in range(0, len(seq), line_length):
            f.write(seq[i: i+line_length] + '\n')
 
def convert_file(input_file, input_type, output_file, output_type):
    """
    Converts files using BioPython AlignIO module.
    """
    with open(input_file, 'r') as f1, open(output_file, 'w') as f2:
        seq = AlignIO.parse(f1, input_type)
        AlignIO.write(seq, f2, output_type)
