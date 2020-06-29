def consensus_variants(vcf_file):
    """
    Inputs:
        vcf_file: a vcf file
    Outputs:
        variants: most common variant that replaces reference in consensus, in the form (position, variant)
    """
    vcf_read = vcf.Reader(open(vcf_file, 'r'))
    variants = []
    previous = (0, 0, 0) #(position, variant, allele frequency)
    previous_previous = (0, 0, 0) #(position, variant, allele frequency)
    for record in vcf_read:
        if record.INFO['AF'] > 0.5: #keeps reference if exactly 0.5
            variants.append((record.POS, record.ALT))
        #2 variants
        if record.POS == previous[0] and record.POS != previous_previous[0]:
            two_ref = (1 - record.POS - previous[2])
            if previous[2] > 0.5:
                previous_previous = previous
                previous = (record.POS, record.ALT, record.INFO['AF'])
                break
            elif two_ref > record.POS and two_ref > previous[0]:
                previous_previous = previous
                previous = (record.POS, record.ALT, record.INFO['AF'])
                break
            elif record.POS > two_ref and record.POS > previous[2]:
                variants.append((record.POS, record.ALT))
            elif previous[2] > two_ref and previous[2] > record.POS:
                variants.append((previous[0], previous[1]))

        #3 variants
        if record.POS == previous[0] and record.POS == previous_previous[0]:
            three_ref = (1 - record.POS - previous[2] - previous_previous[2])
            if previous[2] > 0.5 or previous_previous[2] > 0.5:
                previous_previous = previous
                previous = (record.POS, record.ALT, record.INFO['AF'])
                break
            elif three_ref > record.POS and three_ref > previous[2] and three_ref > previous_previous[2]:
                #ref allele is greatest freq
                previous_previous = previous
                previous = (record.POS, record.ALT, record.INFO['AF'])
                break
            elif record.POS > three_ref and record.POS > previous[2] and record.POS > previous_previous[2]:
                variants.append((record.POS, record.ALT))
            elif previous[2] > two_ref and previous[2] > record.POS and previous[2] > previous_previous[2]:
                variants.append(previous[0], previous[1])
            elif previous_previous[2] > two_ref and previous_previous[2] > record.POS and previous_previous[2] > previous[2]:
                variants.append(previous_previous[0], previous_previous[1])
        previous_previous = previous
        previous = (record.POS, record.ALT, record.INFO['AF'])
    return variants




def consensus(reference, variants):
    """
    Input:
        variants: most common variant that replaces reference in consensus, in the form (position, variant)
        reference: the reference sequence
    Output:
        reference: reference is updated to be the consensus
    """
    for variant in variants:
        reference[(variant[0]) - 1] = variant[1]
    return reference
