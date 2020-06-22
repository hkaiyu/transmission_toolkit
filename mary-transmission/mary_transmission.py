import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import vcf
import subprocess
import os

def lfv_extractor(filename, extract_type):
    """
    Inputs:
    filename: a vcf file
    extract_type: a string that is either 'biallelic' or 'multiallelic'
    Outputs: 
    a list of tuples with low frequency variants in the form (position, variant, frequency)
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
        
    
    

def shared_variants(lfv_list1, lfv_list2):
    """
    Inputs:
        lfv_dict1: a list of tuples containing low frequency variants in the form (position, variant, frequency)
        lfv_dict2: a list of tuples containing low frequency variants in the form (position, variant, frequency)
    Outputs: 
        shared: a list of tuples with all of the shared variants between the two input files in the form (position, variant, frequency1, frequency2)
    """
    #adds filenames as first two elements of shared
    shared = [lfv_list1[0], lfv_list2[0]]
    #get list of shared variants
    for num1 in range(1, len(lfv_list1) - 1):
        for num2 in range(1, len(lfv_list2) - 1):
            #add to list of shared variants
            if lfv_list1[num1][0] == lfv_list2[num2][0] and lfv_list1[num1][1] == lfv_list2[num2][1]:
                shared.append((lfv_list1[num1][0], lfv_list1[num1][1], lfv_list1[num1][2], lfv_list2[num2][2]))
    return shared

def bar_plot_input(donor, recipient):
    """
    Inputs:
        donor: a list of tuples containing low frequency variants in the form (position, variant, frequency)
        recipient: a list of tuples containing low frequency variants in the form (position, variant, frequency)
    Outputs:
        bar_input: a list of tuples with all of the donor variants and their frequencies in the recipient in the form (position, variant, donor frequency, recipient frequency)
    """
    bar_input = [donor[0], recipient[0]]
    for num1 in range(1, len(donor) - 1):
        recipient_frequency = 0
        for num2 in range(1, len(recipient) - 1):
            if donor[num1][0] == recipient[num2][0] and donor[num1][1] == recipient[num2][1]:
                recipient_frequency = recipient[num2][2]
        bar_input.append((donor[num1][0], donor[num1][1], donor[num1][2], recipient_frequency))
    return bar_input
            



def bar_plot(shared_variants):
    """
    Inputs:
    a list of tuples of shared variants between two vcf files in the form (position, variant, frequency1, frequency2)
    Outputs:
    a bar plot containing the shared variants and frequencies of the two input vcf files
    """
    
    filename1 = shared_variants[0]
    filename2 = shared_variants[1]

    positions = []
    frequency1 = []
    frequency2 = []
    allele = []

    for num in range(2, len(shared_variants)):
        positions.append((shared_variants[num][0], shared_variants[num][1]))
        allele.append(shared_variants[num][1])
        frequency1.append(round(shared_variants[num][2], 2))
        frequency2.append(round(shared_variants[num][3], 2))

    x = np.arange(len(positions)) #label location
    width = 0.25 #bar width

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, frequency1, width, label=filename1, color='#7f6d5f')
    rects2 = ax.bar(x + width/2, frequency2, width, label=filename2, color='#557f2d')

    ax.set_ylabel('Allele Frequency')
    ax.set_xlabel('Position & Allele')
    ax.set_title('Low Frequency Variants')
    ax.set_xticks(x)
    ax.set_xticklabels(positions)
    ax.legend()

    plt.xticks(rotation=60, size=9)
    plt.yticks(size=9)
    plt.legend(fontsize=9)

    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height), xy=(rect.get_x() + rect.get_width() / 2, height), xytext=(0, 3), textcoords="offset points", ha='center', va='bottom')
    
    autolabel(rects1)
    autolabel(rects2)

    fig.tight_layout()
    plt.show()
    


file1 = lfv_extractor('COV-20200312-P2-E01-N_S31_bwamem.bam.lowfreq.vcf', 'multiallelic')
file2 = lfv_extractor('COV-20200312-P2-E03-N_S37_bwamem.bam.lowfreq.vcf', 'multiallelic')

#print(file1[2][0])
#print(file2)

prac = (bar_plot_input(file1, file2))

#print(prac)

#bar_plot(prac)


def bottleneck_input(bb_input):
    """
    Inputs:
        bb_input: a list of tuples with all of the donor variants and their frequencies in the recipient in the form (position, variant, donor frequency, recipient frequency)
    Output:
        bottleneck.txt: a txt file in the form of a bb_bottleneck input file
    """
    bb = open("bottleneck.txt", "w+")
    for num in range(2, len(bb_input) - 1):
        bb.write(str(bb_input[num][2]) + "\t" + str(bb_input[num][3]) + "\n")
    bb.close()

bottleneck_input(prac)

(subprocess.run('Rscript Bottleneck_size_estimation_approx.r --file "bottleneck.txt" --plot_bool TRUE --var_calling_threshold 0.03 --Nb_min 1 --Nb_max 200 --confidence_level .95'))


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


   


