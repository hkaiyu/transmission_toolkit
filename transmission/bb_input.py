"""Module that contains functions relating to BB_bottleneck software."""

#Local import
from transmission_toolkit import parsers

def bb_file_writer(donor, recipient, parse_type="biallelic", min_read_depth=0, max_AF=1):
    """
    Writes input file to BB_bottleneck software.
    """
    # Parse data and store outputs
    vcf_data, shared_count = parsers.bb_input_data(donor, recipient)

    # Create shortened filenames
    fname1 = donor.split("_")[0] 
    fname2 = recipient.split("_")[0]

    # Write data to file in BB_Bottleneck input format
    filename = f"new_{fname1}_{fname2}_thred{shared_count}_complete_nofilter_bbn.txt"
    with open(filename, "w") as f:
        for pos in vcf_data:
            for var in vcf_data[pos]:
                f.write(str(vcf_data[pos][var][0])+'\t'+str(vcf_data[pos][var][1])+'\n')


