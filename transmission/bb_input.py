from parsers import parse_biallelic as biallelic, parse_multiallelic as multiallelic

def write_bb_input_files(donor_filename, recipient_filename, parse_type, min_AF = 0.03, max_AF = 0.5):
        """
    Writes two files that will be used as input for the BB Bottleneck software
    Inputs:
    > donor_filename - a string representation of the name of the donor's vcf file
    > recipient_filename - a string representation of the name of the recipient's vcf file
    > min_AF - the minimum frequency at which we consider the allele to be a low frequency 
        variant; if not specified, it is defaulted to be at 0.03
    > max_AF - the maximum frequency at which we consider the allele to be a low frequency 
        variant; if not specified, it is defaulted to be at 0.5
    """
    vcf_data, shared_count = parse_type(donor_filename, recipient_filename)
    fname1 = donor_filename.split("_")[0] 
    fname2 = recipient_filename.split("_")[0] 
    file1 = open("new_%s_%s_thred%i_complete_nofilter_bbn.txt"%(fname1, fname2, shared_count), "w") #rmr to remove 'new'
    #new_COV-20200312-P2-E01-N_COV-20200312-P2-E03-N_thred2_complete_nofilter_bbn
    for pos in vcf_data:
        for var in vcf_data[pos]:
            if vcf_data[pos][var][0] > min_AF:
                file1.write(str(vcf_data[pos][var][0])+'\t'+str(vcf_data[pos][var][2])+'\n')
    file1.close()
