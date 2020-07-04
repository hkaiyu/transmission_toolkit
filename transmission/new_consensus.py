import vcf
from transmission_toolkit import parsers

def consensus_sequence(vcf_file, reference_sequence):
    """
    """
    variants = parsers.extract_lfv(vcf_file)
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
        reference_sequence[pos - 1] = highest_variants[pos]
    
    return reference_sequence
        

