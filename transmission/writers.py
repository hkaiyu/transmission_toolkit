"""Module that contains functions relating to BB_bottleneck software."""

# Standard library imports
import os
import subprocess

#Local import
from transmission_toolkit import parsers

def write_fasta(consensus_dict, output_file, line_length=80):
	"""
	Given a dictionary mapping id's to consensus sequences, returns 
	fasta file with given information.
	"""
	with open(output_file, 'w') as f:
		for id in consensus_dict:
			f.write(f">{id}")
			consensus = consensus_dict[id]
			lines = [consensus[i:i+line_length] for i in range(
                0, len(consensus), line_length)]
			for line in lines:
				f.write(line)
    
def bb_file_writer(donor, recipient, parse_type="biallelic", min_read_depth=0, max_AF=1, var_calling_threshold=0.03):
    """
    Writes input file to BB_bottleneck software.
    """
    # Parse data and store outputs
    vcf_data, shared_count = parsers.bb_input_data(donor, recipient, store_reference=False)[0], parsers.bb_input_data(donor, recipient, store_reference=False)[1]

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