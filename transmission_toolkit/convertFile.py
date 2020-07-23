"""
Module for converting file types using BioPython's SeqIO module
"""
import os
from Bio import AlignIO

# Unfinished function!
def vcf2fasta(vcf_file, reference, name='example', output_dir='', line_length=80):
    """
    Given a directory of VCF files and a reference sequence, writes a fasta for each VCF file.
    """
    # Check if path exists
    if not os.path.exists(vcf_file):
        raise FileNotFoundError(f"File does not exist: {vcf_file}")
    if not os.path.exists(reference):
        raise FileNotFoundError(f"File does not exist: {reference}")

    # If user specifies output directory
    if output_dir:

        # Check if directory already exists
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
    
    # Figure out how to align the sequences first, then write a fasta file
 
def convert_file(input_file, input_type, output_file, output_type):
    """
    Converts files using BioPython AlignIO module.
    """
    with open(input_file, 'r') as f1, open(output_file, 'w') as f2:
        seq = AlignIO.parse(f1, input_type)
        AlignIO.write(seq, f2, output_type)

