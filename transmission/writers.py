"""Module that contains functions relating to BB_bottleneck software."""

# Standard library imports
import os
import subprocess

# Local import
from .consensus import map_consensus
from .parsers import bb_input_data

def write_fasta(vcf_path, reference, output_file='', extension='fasta', output_dir= "", line_length=80):
    """
    Given a dictionary mapping patients to consensus sequences, returns 
    fasta file with given information.
    """
    # Check if path exists
    if not os.path.exists(vcf_path):
        raise IOError(f"Directory does not exist: {vcf_path}")
    if not os.path.exists(reference):
        raise IOError(f"File does not exist: {reference}")

    # If output_file not specified, use default name
    if not output_file:
        output_file = 'example'
    
    filepath = output_file + '.' + extension 

    # If user specifies output directory
    if output_dir:

        # Check if directory already exists
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        
        # Update path to reflect newly added directory
        filepath = os.path.join(output_dir, filepath)

        # If file already exists in this location, remove it to overwrite it
        if os.path.exists(filepath):
            os.remove(filepath)

    consensus_dict = map_consensus(vcf_path, reference)

    # Write file in specified path
    with open(filepath, 'w') as f:
        for fname in consensus_dict:
            f.write(f">{fname}")
            genome = consensus_dict[fname]
            for i in range(0, len(genome), line_length):
                f.write(genome[i: i+line_length])
            
write_fasta('home/kaiyu/Documents/mason_data', 'home/kaiyu/Documents/sequence.fasta')
    
def bb_file_writer(donor, recipient, parse_type="biallelic", min_read_depth=0, max_AF=1, var_calling_threshold=0.03):
    """
    Writes input file to BB_bottleneck software.
    """
    # Parse data and store outputs
    parsed_data = bb_input_data(
        donor,
        recipient,
        parse_type=parse_type,
        min_read_depth=min_read_depth,
        max_AF=max_AF)
        
    vcf_data, shared_count = parsed_data[0], parsed_data[1]

    # Create shortened filenames
    fname1 = donor.split("_")[0] 
    fname2 = recipient.split("_")[0]

    # Write data to file in BB_Bottleneck input format
    filename = f"new_{fname1}_{fname2}_thred{shared_count}_complete_nofilter_bbn.txt"
    with open(filename, "w") as f:
        for pos in vcf_data:
            for var in vcf_data[pos]:
                if vcf_data[pos][var][0] > var_calling_threshold:
                    f.write(str(vcf_data[pos][var][0])+'\t'+str(vcf_data[pos][var][1])+'\n')

