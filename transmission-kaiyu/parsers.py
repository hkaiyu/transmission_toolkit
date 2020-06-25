import vcf

def parse_multiallelic(donor_filename, recipient_filename, site_depth_req, min_AF = 0.03, max_AF = 0.5):
    """
    Extracts low frequency variant data from VCF donor and recipient files and stores data in dictionary.
    Inputs:
    > donor_filename - a string representation of the name of the donor's vcf file
    > recipient_filename - a string representation of the name of the recipient's vcf file
    > site_depth_req - minimum number of reads at a specific site that support a specific variant
    > min_AF - 
    > max_AF - 
    Output:
        Returns a tuple with the first element as a dictionary that stores all of the data needed for BB_Bottleneck software input and
        a donor-recipient allele frequency bar chart, and the second element being an integer representing the number of shared variants
    """
    final_data = [] #(pos, var, donor freq, recip freq)
    reference = (0, 0) #(allele, frequency)
    previous_pos = None
    donor_data = vcf.Reader(open(donor_filename, 'r'))
    for record in donor_data:
        site_depth = record.INFO['DP4'][2] + record.INFO['DP4'][3]
        total_depth = record.INFO['DP4'][0] + record.INFO['DP4'][1]
        freq = site_depth/total_depth
        if previous_pos != record.POS and reference <= 0.5:
            final_data.append(previous_pos, reference[0], reference[1])
        if previous_pos != record.POS:
            reference = (record.REF, (1 - freq))
        if previous_pos == record.POS:
            reference -= freq
        
        #meets minimums to not be error
        if freq >= min_AF and site_depth >= site_depth_req:
            if freq <= max_AF:
                pos, var = record.POS, str(record.ALT[0])
                final_data.append((pos, var, freq))

    

    recipient_data = vcf.Reader(open(recipient_filename, 'r'))
    shared_count = 0
    for record in recipient_data:
        site_depth = record.INFO['DP4'][2] + record.INFO['DP4'][3]
        total_depth = record.INFO['DP4'][0] + record.INFO['DP4'][1]
        freq = site_depth/total_depth
        pos, var = record.POS, str(record.ALT[0])
        for item in final_data:
            if pos == item[0] and var == item[1]:
                shared_count += 1
                item[3] = freq
    
    return (final_data, shared_count)

def parse_biallelic(donor_filename, recipient_filename, site_depth_req, min_AF = 0.03, max_AF = 0.5):
    """
    Extracts low frequency variant data from VCF donor and recipient files and stores data in dictionary.
    Inputs:
    > donor_filename - a string representation of the name of the donor's vcf file
    > recipient_filename - a string representation of the name of the recipient's vcf file
    Output:
        Returns a tuple with the first element as a dictionary that stores all of the data needed for BB_Bottleneck software input and
        a donor-recipient allele frequency bar chart, and the second element being an integer representing the number of shared variants
    """
    final_data = [] #(pos, var, donor freq, recip freq)
    reference = (0, 0)#(allele, frequency)
    previous_pos = None
    donor_data = vcf.Reader(open(donor_filename, 'r'))
    for record in donor_data:
        site_depth = record.INFO['DP4'][2] + record.INFO['DP4'][3]
        total_depth = record.INFO['DP4'][0] + record.INFO['DP4'][1]
        freq = site_depth/total_depth
        if previous_pos != record.POS and reference <= 0.5:
            final_data.append(previous_pos, reference[0], reference[1])
        if previous_pos != record.POS:
            reference = (record.REF, (1 - freq))
        if previous_pos == record.POS:
            reference -= freq
        #meets minimums to not be error
        if freq >= min_AF and site_depth >= site_depth_req:
            if freq <= max_AF:
                pos, var = record.POS, str(record.ALT[0])
                final_data.append((pos, var, freq))

    

    recipient_data = vcf.Reader(open(recipient_filename, 'r'))
    shared_count = 0
    for record in recipient_data:
        site_depth = record.INFO['DP4'][2] + record.INFO['DP4'][3]
        total_depth = record.INFO['DP4'][0] + record.INFO['DP4'][1]
        freq = site_depth/total_depth
        pos, var = record.POS, str(record.ALT[0])
        for item in final_data:
            if pos == item[0] and var == item[1]:
                shared_count += 1
                item[3] = freq
    
    return (final_data, shared_count)
