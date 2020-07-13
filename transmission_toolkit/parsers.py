"""Module containing classes and functions for parsing different types of data"""
# Standard library import
import os

#Local import
from .utils import Data, Biallelic, Multiallelic

#Third party import
import vcf

def _is_valid_lfv(min_read_depth, max_AF, var_reads, total_reads):
    """
    Boolean function called as a helper for determining if a variant fits the
    user-defined parameters for parsing the data.
    """
    # Calculate allele frequency
    freq = var_reads / total_reads

    #If the allele passes restrictions, return True
    if var_reads >= min_read_depth and freq < max_AF:
        return True

    #Otherwise, return False
    return False

def mask_parse(masks):
    """
    Parses a txt file of masked positions
    """
    mask_positions = []
    with open(masks, "r") as mask_file:
            for line in mask_file:
                comma_split = line.split(",")
                for item in comma_split:
                    nums = item.split("-")
                    if len(nums) == 1:
                        mask_positions.append(int(nums[0]))
                    else:
                        for num in range(int(nums[0]), int(nums[1])):
                            mask_positions.append(int(num))
    return mask_positions

def extract_lfv(filepath, min_read_depth=0, max_AF=1, parse_type="biallelic", store_reference=True, masks=None, mask_status='hide'):
    """
    Extracts variant data from VCF and creates a dictionary storing data
    in the form: {position: {variant: [frequency, depth]}}.
    """

    #### Handle Errors #####
    PARSE_TYPES = {"biallelic", "multiallelic"}
    MASK_TYPES = {"hide", "highlight"}
    if min_read_depth < 0 or max_AF > 1 or parse_type not in PARSE_TYPES or str(mask_status) not in MASK_TYPES:
        raise ValueError("Invalid input.")
    #########################

    #Parse mask file if mask file is inputted
    if masks != None:
        mask_positions = mask_parse(masks)

    lfv_data = Biallelic() if parse_type == "biallelic" else Multiallelic()
    data = vcf.Reader(open(filepath, 'r'))
    ref_data = {}

    for record in data:

        var_depth = record.INFO['DP4'][2] + record.INFO['DP4'][3] # Num of variant reads
        raw_depth = record.INFO['DP'] # Num of total reads
        pos, var, freq = record.POS, str(record.ALT[0]), float(var_depth / raw_depth) 

        #doesn't include masked positions based on user settings
        if masks!= None and pos in mask_positions and mask_status == 'hide':
            continue

        # If variant passes restrictions, store data
        if _is_valid_lfv(min_read_depth, max_AF, var_depth, raw_depth):
            lfv_data.store(pos, var, freq, var_depth)

        if store_reference and not pos in ref_data:

            # Get reference data
            ref_depth = record.INFO['DP4'][0] + record.INFO['DP4'][1]
            ref = str(record.REF[0])
            
            # If ref allele passes restrictions, store the data
            if _is_valid_lfv(min_read_depth, max_AF, ref_depth, raw_depth):
                ref_data[pos] = {ref: [(ref_depth / raw_depth), ref_depth]}

    # After parsing is complete, make object into a dictionary
    lfv_data = dict(lfv_data)

    # If we collected reference data, update lfv_data
    if ref_data:
        for pos in ref_data:
            if not pos in lfv_data:
                lfv_data[pos] = ref_data[pos]
            else:
                lfv_data[pos].update(ref_data[pos])

    return lfv_data

def bb_input_data(donor, recip, min_read_depth=0, max_AF=1, parse_type="biallelic", store_reference=True, weighted=False):
    """
    Stores info from parsing VCF files to dictionary.
    """
    donor_data = extract_lfv(donor, min_read_depth=min_read_depth, max_AF=max_AF, parse_type=parse_type, store_reference=store_reference)
    recip_data = extract_lfv(recip, min_read_depth=min_read_depth, max_AF=max_AF, parse_type=parse_type, store_reference=store_reference)

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


    return (bb_data, shared_count)
