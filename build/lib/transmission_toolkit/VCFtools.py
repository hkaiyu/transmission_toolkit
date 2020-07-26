"""Module for extracting data from VCF files"""
import os
import glob
import vcf
from transmission_toolkit.utils import get_seq

def _is_lfv(min_read_depth, max_AF, var_reads, total_reads):
    freq = var_reads / total_reads
    if var_reads >= min_read_depth and freq < max_AF:
        return True
    return False

def _biallelic_store(var_dict, pos, var, freq, var_depth):
    if pos in var_dict:
        variant = list(var_dict[pos].values())
        old_freq = variant[0][0]
        if old_freq < freq:
            var_dict[pos] = {var: [freq, var_depth]}
    else:
        var_dict[pos] = {var: [freq, var_depth]}

def _multiallelic_store(var_dict, pos, var, freq, var_depth):
        if pos in var_dict:
            var_dict[pos][var] = [freq, var_depth]
        else:
            var_dict[pos] = {var: [freq, var_depth]}

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

def extract_lfv(
    vcf_file, 
    min_read_depth=1, 
    max_AF=1, 
    parse_type='biallelic', 
    store_ref=True, 
    masks=None, 
    mask_status='hide'
    ):
    """
    Extracts variant data from VCF and creates a dictionary storing data
    in the form: {position: {variant: [frequency, depth]}}.
    """

    #### Handle Errors #####
    PARSE_TYPES = {"biallelic", "multiallelic"}
    MASK_TYPES = {"hide", "highlight"}
    if min_read_depth < 0 or max_AF > 1 or parse_type not in PARSE_TYPES or mask_status not in MASK_TYPES:
        raise ValueError("Invalid input.")
    #########################

    #Parse mask file if mask file is inputted
    if masks != None:
        mask_positions = mask_parse(masks)

    lfv_data, ref_data = {}, {}
    data = vcf.Reader(open(vcf_file, 'r'))

    for record in data:

        var_depth = record.INFO['DP4'][2] + record.INFO['DP4'][3] # Num of variant reads
        raw_depth = record.INFO['DP'] # Num of total reads
        pos, var, freq = record.POS, str(record.ALT[0]), float(var_depth / raw_depth) 

        #doesn't include masked positions based on user settings
        if masks!= None and pos in mask_positions and mask_status == 'hide':
            continue

        # If variant passes restrictions, store data
        if _is_lfv(min_read_depth, max_AF, var_depth, raw_depth):
            if parse_type == 'biallelic':
                _biallelic_store(lfv_data, pos, var, freq, var_depth)
            else:
                _multiallelic_store(lfv_data, pos, var, freq, var_depth)

        if store_ref and not pos in ref_data:

            # Get reference data
            ref_depth = record.INFO['DP4'][0] + record.INFO['DP4'][1]
            ref = str(record.REF[0])
            
            # If ref allele passes restrictions, store the data
            if _is_lfv(min_read_depth, max_AF, ref_depth, raw_depth):
                ref_data[pos] = {ref: [(ref_depth / raw_depth), ref_depth]}

    # If we collected reference data, update lfv_data
    if ref_data:
        for pos in ref_data:
            if not pos in lfv_data:
                lfv_data[pos] = ref_data[pos]
            else:
                lfv_data[pos].update(ref_data[pos])

    return lfv_data

def build_majority_consensus(
    vcf_file, 
    reference,  
    store_ref=True, 
    masks=None, 
    mask_status='hide'
    ):
    """
    With a reference file specified, builds consensus sequence.
    """
    
    consensus = list(get_seq(reference))
    variants = extract_lfv(
        vcf_file, 
        min_read_depth=1, 
        max_AF=1, 
        parse_type='biallelic', 
        store_ref=store_ref, 
        masks=None, 
        mask_status='hide'
    )

    highest_variants = {}

    for pos in variants:
        highest_freq = 0
        kept_var = None
        for var in variants[pos]:
            if variants[pos][var][0] > highest_freq:
                highest_freq = variants[pos][var][0]
                kept_var = var
        if kept_var:
            highest_variants[pos] = str(kept_var)

    for pos in highest_variants:
        consensus[pos - 1] = highest_variants[pos]

    return ''.join(consensus)

def build_minor_consensus(
    vcf_file, 
    reference, 
    min_read_depth=1, 
    max_AF=1, 
    store_ref=True, 
    masks=None, 
    mask_status='hide'
    ):
    """
    Builds consensus with variants that have second highest frequency.
    """
    consensus = list(get_seq(reference))
    data = extract_lfv(
        vcf_file, 
        min_read_depth=min_read_depth, 
        max_AF=max_AF, 
        parse_type='multiallelic', 
        store_ref=store_ref, 
        masks=None, 
        mask_status='hide'
    )

    variants = {}

    for pos in data:
        highest_freq = 0
        max_freq = max([data[pos][var][0] for var in data[pos]])
        for var in data[pos]:
            freq = data[pos][var][0]
            if freq > highest_freq and freq != max_freq:
                highest_freq = data[pos][var][0]
                kept_var = var
        variants[pos] = str(kept_var)

    for pos in variants:
        consensus[pos - 1] = variants[pos]

    return ''.join(consensus)


#data = VCFtools('example_data/test1.vcf')
#print(data.extract_lfv(store_reference=False))
#print(data.build_consensus('example_data/testref.fasta'))