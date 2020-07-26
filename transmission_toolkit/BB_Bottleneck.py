"""Module containing classes and functions for parsing different types of data"""
# Standard library import
import os

#Local import
from transmission_toolkit.VCFtools import extract_lfv, build_majority_consensus

def bb_input_data(donor, recip, min_read_depth=0, max_AF=1, parse_type="biallelic",  store_ref=True, weighted=False, masks=None, mask_status="hide"):
    """
    Stores info from parsing VCF files to dictionary.
    """
    donor_data, donor_masks = extract_lfv(
        donor, 
        min_read_depth=min_read_depth, 
        max_AF=max_AF, 
        parse_type=parse_type, 
        store_ref=store_ref, 
        masks=masks, 
        mask_status=mask_status
    )
    recip_data, recip_masks = extract_lfv(
        recip, 
        min_read_depth=min_read_depth, 
        max_AF=max_AF, 
        parse_type=parse_type, 
        store_ref=store_ref, 
        masks=masks, 
        mask_status=mask_status
    )

    shared_count = 0

    # Stored as {pos: {var: [donor freq., recip. freq]}} bc Maria had two bb input files 
    # and one required all this info, might change later tho
    bb_data = {} 

    # Iterate through each variant at each position
    for pos in donor_data: 

        bb_data[pos] = {}

        for var in donor_data[pos]:

            # Store donor data
            donor_freq = donor_data[pos][var][0]
            bb_data[pos] = {var: [donor_freq, 0.0]}
            if weighted == True:
                donor_depth = donor_data[pos][var][1]
                bb_data[pos] = {var: [donor_freq, 0.0, donor_depth, 0.0]}

            # If recipient has same variant at same location, store it
            if pos in recip_data and var in recip_data[pos]:
                recip_freq = recip_data[pos][var][0]
                bb_data[pos][var][1] = recip_freq
                shared_count += 1
                if weighted == True:
                    recip_depth = recip_data[pos][var][1]
                    bb_data[pos][var][3] = recip_depth


    return (bb_data, shared_count, donor_masks, recip_masks)
    
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