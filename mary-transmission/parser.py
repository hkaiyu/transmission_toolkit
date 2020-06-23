import vcf

def lfv_extractor(filename, extract_type):
    """
    Inputs:
    filename: a vcf file
    extract_type: a string that is either 'biallelic' or 'multiallelic'
    Outputs: 
    a list of tuples with low frequency variants in the form (position, variant, frequency)
    *filename is first element of output tuple
    *If a position has multiple variants, it will be reflected in multiple tuples
    """
    if extract_type != 'multiallelic' and extract_type != 'biallelic':
        return('Please enter valid extract type.')
    lfv_list = [filename]
    vcf_read = vcf.Reader(open(filename, 'r'))
    previous_previous = (None, None, None)
    previous = (None, None, None)
    for record in vcf_read:
        af = None
        variant = None
        #removes multiallelic variants at a single position
        if previous == record.POS and extract_type == 'biallelic':
            lfv_list.pop()
            previous = record.POS
            break
        #handles case of 2 variants
        if previous[0] == record.POS and previous_previous[0] != record.POS and extract_type == 'multiallelic':
            if previous[1] == record.REF:
                new_af = previous[2] - record.INFO['AF']
                previous = (previous[0], previous[1], new_af)
                lfv_list.pop()
                lfv_list.append(previous)
                lfv_list.append((record.POS, record.ALT, record.INFO['AF']))
                previous_previous = previous
                previous = (record.POS, record.ALT, record.INFO['AF'])
            elif record.INFO['AF'] > 0.5:
                af = 1 - record.INFO['AF'] - previous[2]
                lfv_list.append((record.POS, record.REF, af))
                previous_previous = previous
                previous = (record.POS, record.REF, af)
            else: 
                lfv_list.append((record.POS, record.ALT, record.INFO['AF']))
                previous_previous = previous
                previous = (record.POS, record.ALT, record.INFO['AF'])
            break
        if previous[0] == record.POS and previous_previous[0] == record.POS and extract_type == 'multiallelic':
            if previous_previous[1] == record.REF:
                previous_previous[2] = previous_previous[2] - record.INFO['AF']
                lfv_list.pop()
                lfv_list.pop()
                lfv_list.append(previous_previous)
                lfv_list.append(previous)
                lfv_list.append((record.POS, record.ALT, record.INFO['AF']))
                previous_previous = previous
                previous = (record.POS, record.ALT, record.INFO['AF'])
            elif previous[1] == record.REF:
                new_af = previous[2] - record.INFO['AF']
                previous = (previous[0], previous[1], new_af)
                lfv_list.pop()
                lfv_list.append(previous)
                lfv_list.append((record.POS, record.ALT, record.INFO['AF']))
                previous_previous = previous
                previous = (record.POS, record.REF, record.INFO['AF'])
            elif record.INFO['AF'] > 0.5:
                af = 1 - record.INFO['AF'] - previous[2] - previous_previous[2]
                lfv_list.append((record.POS, record.REF, af))
                previous_previous = previous
                previous = (record.POS, record.REF, af)
            else:
                lfv_list.append(record.POS, record.ALT, record.INFO['AF'])
                previous_previous = previous
                previous = (record.POS, record.ALT, record.INFO['AF'])
            break
        #handles cases where the current and previous position are different
        elif record.INFO['AF'] > 0.5:
            af = 1 - record.INFO['AF']
            variant = record.REF
        elif record.INFO['AF'] < 0.5:
            af = record.INFO['AF']
            variant = record.ALT
        lfv_list.append((record.POS, variant, af))
        previous_previous= previous
        previous = (record.POS, variant, af)
    return lfv_list