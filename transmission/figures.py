import matplotlib.pyplot as plt
import datetime
import numpy as np
from pathlib import Path
import os
from transmission_toolkit import parsers

def make_standard_bar_plot(donor_filename, recipient_filename, min_read_depth=0, max_AF=1, parse_type="biallelic", store_reference=True, plot_type='standard', masks='default_mask', mask_type='hide'):
    """
    Saves a barplot figure showing allele frequencies vs. position.
    Inputs:
    > donor_filename - a string representation of the name of the donor's vcf file
    > recipient_filename - a string representation of the name of the recipient's vcf file
    > min_AF - the minimum frequency at which we consider the allele to be a low frequency 
        variant; if not specified, it is defaulted to be at 0.03
    > max_AF - the maximum frequency at which we consider the allele to be a low frequency 
        variant; if not specified, it is defaulted to be at 0.5
    """
    PARSE_TYPES = {"biallelic", "multiallelic"}
    PLOT_TYPES = {"standard", "weighted"}
    MASK_TYPES = {"hide", "highlight"}
    if parse_type not in PARSE_TYPES or plot_type not in PLOT_TYPES or mask_type not in MASK_TYPES:
        raise ValueError("Invalid input.")


    vcf_data = parsers.bb_input_data(donor_filename, recipient_filename, min_read_depth=min_read_depth, max_AF=max_AF, parse_type=parse_type, store_reference=store_reference, weighted=True, masks=masks, mask_type=mask_type)[0]
    positions, donor_freqs, recipient_freqs, donor_depths, recipient_depths = [], [], [], [], []

    donor_filename = donor_filename.split("_")
    recipient_filename = recipient_filename.split("_")
    donor_label = str('Donor: ' + donor_filename[0])
    recipient_label = str('Recipient: ' + recipient_filename[0])

    for pos in sorted(vcf_data):
        for nuc in vcf_data[pos]:
            positions.append((pos, nuc))
            donor_freqs.append(round(vcf_data[pos][nuc][0], 2))
            recipient_freqs.append(round(vcf_data[pos][nuc][1], 2))
            donor_depths.append((vcf_data[pos][nuc][2])/50)
            recipient_depths.append((vcf_data[pos][nuc][3])/50)

    x = np.arange(len(positions))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()

    if plot_type == 'standard':
        rects1 = ax.bar(x - width/2, donor_freqs,  color='#7f6d5f', width=width, \
            edgecolor='white', label=donor_label)
        rects2 = ax.bar(x + width/2, recipient_freqs, color='#557f2d', width=width, \
            edgecolor='white', label=recipient_label)

    if plot_type == 'weighted':
        rects1 = ax.bar(x - width/2, donor_freqs,  color='#7f6d5f', width=donor_depths, \
            edgecolor='white', label=donor_label)
        rects2 = ax.bar(x + width/2, recipient_freqs, color='#557f2d', width=recipient_depths, \
            edgecolor='white', label=recipient_label)
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Allele Frequencies')
    ax.set_xticks(x)
    ax.set_xticklabels(positions)
    ax.tick_params(axis = 'x', rotation = 50)
    ax.legend()

    
    maxn = len(positions)
    fig.tight_layout() 
    plt.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2, borderaxespad=0, frameon=False, fontsize=(maxn/1.6))
    plt.xlim(-0.5,maxn-0.5)
    plt.tight_layout()
    plt.savefig('%s_%s_allele_freq.png'%(donor_filename, recipient_filename), dpi=300, bbox_inches='tight')
    plt.show()


#make_standard_bar_plot('COV-20200312-P2-E01-N_S31_bwamem.bam.lowfreq.vcf', 'COV-20200312-P2-E03-N_S37_bwamem.bam.lowfreq.vcf', plot_type='weighted')

def masked_shared_variants(vcf_path, masks, min_read_depth=0, max_AF=1, parse_type="biallelic", store_reference=True):
    """
    {number of shared variants: [complete, masked]}
    """
    vcf_folder = Path(vcf_path)
    shared_variants = {}
    for donor_filename in os.listdir(vcf_folder):
        for recipient_filename in os.listdir(vcf_folder):
            if donor_filename != recipient_filename:
                masked_shared_count = parsers.bb_input_data(donor_filename, recipient_filename, store_reference=True, masks=masks, mask_type='hide')[1]
                complete_shared_count = parsers.bb_input_data(donor_filename, recipient_filename, store_reference=True, masks=masks, mask_type='highlight')[1]


                if complete_shared_count in shared_variants:
                    shared_variants[complete_shared_count][0] += 1
                elif complete_shared_count not in shared_variants:
                    shared_variants[complete_shared_count] = [1, 0]

                if masked_shared_count in shared_variants:
                    shared_variants[masked_shared_count][1] += 1
                elif masked_shared_count not in shared_variants:
                    shared_variants[masked_shared_count] = [0, 1]

    shared_variants = sorted(shared_variants.items())
    sv, masked_pairs, complete_pairs = [], [], []
    for num in shared_variants:
        sv.append(num[0])
        complete_pairs.append(num[1][0])
        masked_pairs.append(num[1][1])


    masked_label = "Masked Genome"
    complete_label = "Complete Genome"

    x = np.arange(len(sv))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()


    rects1 = ax.bar(x - width/2, masked_pairs,  color='red', width=width, \
        edgecolor='white', label=masked_label)
    rects2 = ax.bar(x + width/2, complete_pairs, color='blue', width=width, \
        edgecolor='white', label=complete_label)
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Number of Pairs')
    ax.set_xlabel('Number of Shared Variants')
    ax.set_xticks(x)
    ax.set_xticklabels(sv)
    ax.tick_params(axis = 'x')
    ax.legend()

    
    maxn = len(sv)
    fig.tight_layout() 
    plt.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2, borderaxespad=0, frameon=False, fontsize=(maxn/1.6))
    plt.xlim(-0.5,maxn-0.5)
    plt.xticks(fontsize=8)
    #plt.tight_layout()
    #plt.savefig('%s_%s_allele_freq.png'%(donor_filename, recipient_filename), dpi=300, bbox_inches='tight')
    plt.show()
