import vcf

def parse_multiallelic(donor_filename, recipient_filename, min_AF = 0.03, max_AF = 0.5):
    """
    Extracts low frequency variant data from VCF donor and recipient files and stores data in dictionary.

    Inputs:
    > donor_filename - a string representation of the name of the donor's vcf file
    > recipient_filename - a string representation of the name of the recipient's vcf file

    Output:
        Returns a tuple with the first element as a dictionary that stores all of the data needed for BB_Bottleneck software input and
        a donor-recipient allele frequency bar chart, and the second element being an integer representing the number of shared variants
    """
    final_data = {} #{pos: {var: [donor freq, recip freq]}}
    reference_alleles = {} #{pos: reference}

    donor_data = vcf.Reader(open(donor_filename, 'r'))
    for record in donor_data:
        freq = record.INFO['AF']
        if freq >= min_AF and freq <= max_AF:
            pos, var = record.POS, str(record.ALT[0])
            if pos in final_data:
                final_data[pos][var] = [freq, 0.0]
            else:
                final_data[pos] = {var: [freq, 0.0]}
            reference_alleles[pos] = str(record.REF[0])

    for pos in donor_data:
        total = 0
        for freq in final_data[pos].values():
            total += freq[0]
        final_data[pos][reference_alleles[pos]] = 1 - total

    recipient_data = vcf.Reader(open(recipient_filename, 'r'))
    shared_count = 0
    for record in recipient_data:
        pos, var = record.POS, str(record.ALT[0])
        if pos in final_data and var in final_data[pos]:
            shared_count += 1
            final_data[pos][var][1] = record.INFO['AF']
    
    return (final_data, shared_count)

def parse_biallelic(donor_filename, recipient_filename, min_AF = 0.03, max_AF = 0.5):
    """
    Extracts low frequency variant data from VCF donor and recipient files and stores data in dictionary.

    Inputs:
    > donor_filename - a string representation of the name of the donor's vcf file
    > recipient_filename - a string representation of the name of the recipient's vcf file

    Output:
        Returns a tuple with the first element as a dictionary that stores all of the data needed for BB_Bottleneck software input and
        a donor-recipient allele frequency bar chart, and the second element being an integer representing the number of shared variants
    """
    final_data = {}

    donor_data = vcf.Reader(open(donor_filename, 'r'))
    for record in donor_data:
        freq = record.INFO['AF']
        if freq >= min_AF and freq <= max_AF:
            pos, var = record.POS, str(record.ALT[0]) 
            if pos in final_data:
                old_var = list(final_data[pos].keys())[0] #only one item per position bc biallelic
                if freq > final_data[pos][old_var][0]: #if higher freq. lfv found at this site, change old lfv
                    del final_data[pos][old_var]
                    final_data[pos][var] = [freq, 0.0]
            else:
                final_data[pos] = {var: [freq, 0.0]}
                
    recipient_data = vcf.Reader(open(recipient_filename, 'r'))
    shared_count = 0
    for record in recipient_data:
        pos, var = record.POS, str(record.ALT[0])
        if pos in final_data and var in final_data[pos]:
            shared_count += 1
            final_data[pos][var][1] = record.INFO['AF']
    
    return (final_data, shared_count)




