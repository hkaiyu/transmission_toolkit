import vcf

def parse_multiallelic(donor_filename, recipient_filename, allele_depth_req = 1, min_AF = 0.03, max_AF = 0.5):
    """
    Extracts low frequency variant data from VCF donor and recipient files and stores data in dictionary.
    Inputs:
    > donor_filename - a string representation of the name of the donor's vcf file
    > recipient_filename - a string representation of the name of the recipient's vcf file
    > allele_depth_req - minimum number of reads at a specific allele that support a specific variant
    > min_AF - 
    > max_AF - 
    Output:
        Returns a tuple with the first element as a dictionary that stores all of the data needed for BB_Bottleneck software input and
        a donor-recipient allele frequency bar chart, and the second element being an integer representing the number of shared variants
    """
    final_data = {} #{pos: {var: [donor freq, donor depth, recip freq, recip depth]}}
    prev_reference = [None, None, None] #(allele, frequency, depth)
    previous_pos = None
    donor_data = vcf.Reader(open(donor_filename, 'r'))
    for record in donor_data:
        pos, var= record.POS, str(record.ALT[0])
        allele_depth = record.INFO['DP4'][2] + record.INFO['DP4'][3]
        total_depth = record.INFO['DP4'][0] + record.INFO['DP4'][1] + record.INFO['DP4'][2] + record.INFO['DP4'][3]
        freq = float(allele_depth/total_depth) 

        #meets minimums to not be error
        if freq >= min_AF and allele_depth >= allele_depth_req:
            if freq <= max_AF:
                if pos in final_data:
                    final_data[pos].update({var: [freq, allele_depth, 0, 0]})
                else:
                    final_data[pos] = {var: [freq, allele_depth, 0, 0]}

        #handles reference
        if previous_pos != None and previous_pos != pos and prev_reference[1] <= 0.5:
            if previous_pos in final_data:
                final_data[previous_pos].update({prev_reference[0]: [prev_reference[1], prev_reference[2], 0, 0]})
            else:
                final_data[previous_pos] = {prev_reference[0]: [prev_reference[1], prev_reference[2], 0, 0]}
        if previous_pos != pos:
            prev_reference = [record.REF, float(1 - freq), float(total_depth - allele_depth)]
        if previous_pos == pos:
            prev_reference[1] -= freq
            prev_reference[2] -= allele_depth
        previous_pos = pos
        
    
    #adds final reference if necessary
    if (previous_pos not in final_data or prev_reference[0] not in final_data[previous_pos]) and prev_reference[1] <= 0.5:
        if previous_pos not in final_data:
            final_data[previous_pos] = {prev_reference[0]: [prev_reference[1], prev_reference[2], 0, 0]}
        else:
            final_data[previous_pos].update({prev_reference[0]: [prev_reference[1], prev_reference[2], 0, 0]})
    
    #keeps highest variant at a position (biallelic only)
    for position in final_data:
        largest_var_freq = 0
        if len(final_data[position]) > 1:
            for variant in final_data[position]:
                    if final_data[position][variant][0] > largest_var_freq:
                        largest_var_freq = final_data[position][variant][0]
                        kept_variant = variant
                        kept_depth = final_data[position][variant][1]
            final_data[position] = {kept_variant: [largest_var_freq, kept_depth , 0, 0]}

    recipient_data = vcf.Reader(open(recipient_filename, 'r'))
    shared_count = 0
    for record in recipient_data:
        allele_depth = record.INFO['DP4'][2] + record.INFO['DP4'][3]
        total_depth = record.INFO['DP4'][0] + record.INFO['DP4'][1]
        freq = allele_depth/total_depth
        pos, var = record.POS, str(record.ALT[0])
        for position in final_data:
            if pos == position:
                for variant in final_data[position]:
                    if var == variant:
                        shared_count += 1
                        final_data[position][var][2] = freq
                        final_data[position][var][3] = allele_depth
    
    return (final_data, shared_count)

def parse_biallelic(donor_filename, recipient_filename, allele_depth_req = 1, min_AF = 0.03, max_AF = 0.5):
    """
    Extracts low frequency variant data from VCF donor and recipient files and stores data in dictionary.
    Inputs:
    > donor_filename - a string representation of the name of the donor's vcf file
    > recipient_filename - a string representation of the name of the recipient's vcf file
    Output:
        Returns a tuple with the first element as a dictionary that stores all of the data needed for BB_Bottleneck software input and
        a donor-recipient allele frequency bar chart, and the second element being an integer representing the number of shared variants
    """
    final_data = {} #{pos: {var: [donor freq, donor depth, recip freq, recip depth]}}
    prev_reference = [None, None, None] #(allele, frequency, depth)
    previous_pos = None
    donor_data = vcf.Reader(open(donor_filename, 'r'))
    for record in donor_data:
        pos, var= record.POS, str(record.ALT[0])
        allele_depth = record.INFO['DP4'][2] + record.INFO['DP4'][3]
        total_depth = record.INFO['DP4'][0] + record.INFO['DP4'][1] + record.INFO['DP4'][2] + record.INFO['DP4'][3]
        freq = float(allele_depth/total_depth) 

        #meets minimums to not be error
        if freq >= min_AF and allele_depth >= allele_depth_req:
            if freq <= max_AF:
                if pos in final_data:
                    final_data[pos].update({var: [freq, allele_depth, 0, 0]})
                else:
                    final_data[pos] = {var: [freq, allele_depth, 0, 0]}

        #handles reference
        if previous_pos != None and previous_pos != pos and prev_reference[1] <= 0.5:
            if previous_pos in final_data:
                final_data[previous_pos].update({prev_reference[0]: [prev_reference[1], prev_reference[2], 0, 0]})
            else:
                final_data[previous_pos] = {prev_reference[0]: [prev_reference[1], prev_reference[2], 0, 0]}
        if previous_pos != pos:
            prev_reference = [record.REF, float(1 - freq), float(total_depth - allele_depth)]
        if previous_pos == pos:
            prev_reference[1] -= freq
            prev_reference[2] -= allele_depth
        previous_pos = pos
        
    
    #adds final reference if necessary
    if (previous_pos not in final_data or prev_reference[0] not in final_data[previous_pos]) and prev_reference[1] <= 0.5:
        if previous_pos not in final_data:
            final_data[previous_pos] = {prev_reference[0]: [prev_reference[1], prev_reference[2], 0, 0]}
        else:
            final_data[previous_pos].update({prev_reference[0]: [prev_reference[1], prev_reference[2], 0, 0]})
    
    #keeps highest variant at a position (biallelic only)
    for position in final_data:
        largest_var_freq = 0
        if len(final_data[position]) > 1:
            for variant in final_data[position]:
                    if final_data[position][variant][0] > largest_var_freq:
                        largest_var_freq = final_data[position][variant][0]
                        kept_variant = variant
                        kept_depth = final_data[position][variant][1]
            final_data[position] = {kept_variant: [largest_var_freq, kept_depth , 0, 0]}

    recipient_data = vcf.Reader(open(recipient_filename, 'r'))
    shared_count = 0
    for record in recipient_data:
        allele_depth = record.INFO['DP4'][2] + record.INFO['DP4'][3]
        total_depth = record.INFO['DP4'][0] + record.INFO['DP4'][1]
        freq = allele_depth/total_depth
        pos, var = record.POS, str(record.ALT[0])
        for position in final_data:
            if pos == position:
                for variant in final_data[position]:
                    if var == variant:
                        shared_count += 1
                        final_data[position][var][2] = freq
                        final_data[position][var][3] = allele_depth
    
    return (final_data, shared_count)
