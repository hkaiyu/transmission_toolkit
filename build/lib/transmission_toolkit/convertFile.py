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
    min_read_depth=1,
    max_AF=1,
    parse_type='biallelic',
    store_ref=True,
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
            store_ref=store_ref, 
            masks=None, 
            mask_status='hide'
        )
    elif consensus_type == 'minor':
        seq = build_minor_consensus(
            vcf_path, 
            reference, 
            min_read_depth=1, 
            max_AF=1, 
            store_ref=True, 
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

def root_newick(newick, output_file):
    """
    Writes a newick that treats the reference sequence as a root node rather than child node.
    """
    with open(newick, 'r') as f1, open(output_file, 'w') as f2:
        lines = [line.strip() for line in f1.readlines()]
        parenthesis_count = 0
        start, stop = None, None
        for i in range(len(lines[-1])-1, -1, -1):
            if lines[-1][i] == ')':
                parenthesis_count += 1
            if lines[-1][i] == ',':
                for j in range(i, len(lines[-1])):
                    if lines[-1][j] == ':':
                        start, stop = i, j
                        break
                break
        if start and stop:
            lines[-1] = lines[-1][:start] + ')' * parenthesis_count + lines[-1][start+1:stop] + ';'
        else:
            raise ValueError("File in unexpected format.")
        
        for line in lines:
            f2.write(line)

#root_newick('denmark_wgs_47/parsnp/parsnp.tree', 'demo_newick')   

