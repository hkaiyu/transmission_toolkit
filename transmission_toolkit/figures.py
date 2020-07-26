"""Module for generating bottleneck figures"""

import matplotlib.pyplot as plt
import datetime
import numpy as np
from pathlib import Path
import os
from transmission_toolkit.BB_Bottleneck import bb_input_data
from transmission_toolkit.VCFtools import mask_parse


def bar_plots(vcf_path, masks=None, mask_status='hide', min_read_depth=10, max_AF=1, parse_type='biallelic'):
    all_pairs = []
    vcf_folder = Path(vcf_path)
    for donor_filename in os.listdir(vcf_folder):
        for recipient_filename in os.listdir(vcf_folder):
            if donor_filename != recipient_filename:

def make_standard_bar_plot(
    donor_filename, 
    recipient_filename, 
    min_read_depth=0, 
    max_AF=1, 
    parse_type="biallelic", 
    store_ref=True, 
    plot_type='standard', 
    masks='default_mask', 
    mask_status='hide'
    ):
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
    mask_statusS = {"hide", "highlight"}
    if parse_type not in PARSE_TYPES or plot_type not in PLOT_TYPES or mask_status not in mask_statusS:
        raise ValueError("Invalid input.")


    vcf_data = bb_input_data(
        donor_filename, 
        recipient_filename, 
        min_read_depth=min_read_depth, 
        max_AF=max_AF, 
        parse_type=parse_type, 
        store_ref=store_ref, 
        weighted=True, 
        masks=masks, 
        mask_status=mask_status
    )[0]
    
    positions, donor_freqs, recipient_freqs, donor_depths, recipient_depths = [], [], [], [], []

    donor_filename = donor_filename.split("_")
    recipient_filename = recipient_filename.split("_")
    donor_label = str('Donor: ' + donor_filename[0])
    recipient_label = str('Recipient: ' + recipient_filename[0])
    fig_save = str(donor_filename[0] + recipient_filename[0])

    max_read = 0
    min_read = float('inf')

    for pos in sorted(vcf_data):
        for nuc in vcf_data[pos]:
            positions.append((pos, nuc))
            donor_freqs.append(round(vcf_data[pos][nuc][0], 2))
            recipient_freqs.append(round(vcf_data[pos][nuc][1], 2))
            donor_depths.append((vcf_data[pos][nuc][2])/50)
            recipient_depths.append((vcf_data[pos][nuc][3])/50)
            if vcf_data[pos][nuc][2] > max_read:
                vcf_data[pos][nuc][2] = max_read
            if vcf_data[pos][nuc][3] > max_read:
                vcf_data[pos][nuc][3] = max_read
            if vcf_data[pos][nuc][2] < min_read:
                vcf_data[pos][nuc][2] = min_read
            if vcf_data[pos][nuc][3] < min_read:
                vcf_data[pos][nuc][3] = min_read


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

#masked_shared_variants('mason_data/', 'default_mask.txt', max_AF=0.5, min_read_depth=10)
#make_standard_bar_plot('COV-20200312-P2-E01-N_S31_bwamem.bam.lowfreq.vcf', 'COV-20200312-P2-E03-N_S37_bwamem.bam.lowfreq.vcf', plot_type='weighted')
masked = all_pairs_parse('mason_data/', 'default_mask.txt', max_AF=0.5, min_read_depth=10)
complete = all_pairs_parse('mason_data/', max_AF=0.5, min_read_depth=10)
shared_variants = sv_count(complete, masked)


def masked_shared_variants(sv_count):
    """
    {number of shared variants: [complete, masked]}
    chart c
    """

    shared_variants = sorted(sv_count.items())
    sv, masked_pairs, complete_pairs, cell_text, columns, masked, complete= [], [], [], [], [], [], []
    for num in shared_variants:
        sv.append(num[0])
        complete_pairs.append(num[1][0])
        complete.append(num[1][0])
        masked_pairs.append(num[1][1])
        masked.append(num[1][1])
    
    cell_text.append(masked)
    cell_text.append(complete)

    columns = tuple(sv)
    rows = ['Complete Genome', 'Masked Genome']
    colors = ['blue', 'red']


    masked_label = "Masked Genome"
    complete_label = "Complete Genome"

    x = np.arange(len(sv))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()


    rects1 = ax.bar(x - width/2, masked_pairs,  color='red', width=width, \
        edgecolor='white', label=masked_label)
    rects2 = ax.bar(x + width/2, complete_pairs, color='blue', width=width, \
        edgecolor='white', label=complete_label)

    plt.table(cellText=cell_text, rowLabels=rows, rowColours=colors, colLabels=columns)
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Number of Pairs')
    ax.set_xlabel('Number of Shared Variants')
    #ax.set_xticks(x)
    #ax.set_xticklabels(sv)
    #ax.tick_params(axis = 'x')
    ax.legend()

    
    maxn = len(sv)
    fig.tight_layout()
    #plt.figure(figsize=(10,8)) 
    plt.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2, borderaxespad=0, frameon=False, fontsize=(12))
    plt.xlim(-0.5,maxn-0.5)
    plt.xticks([])
    plt.xlabel("Shared Variants", labelpad=50)
    #plt.tight_layout()
    #plt.savefig('%s_%s_allele_freq.png'%(donor_filename, recipient_filename), dpi=300, bbox_inches='tight')
    plt.show()


#masked_shared_variants('mason_data/', 'default_mask.txt', max_AF=0.5)

def shared_positions(position_count, mask_file):
    """
    """
    mask_positions = mask_parse(mask_file)

    sorted_pos = sorted(position_count.items(), key=lambda x: x[1], reverse=True)
    position_data, frequency_data, colors = [], [], []
    sorted_pos = sorted_pos[0:shown_variants]
    print(sorted_pos)
    for pos in sorted_pos:
        position_data.append(pos[0])
        if pos[0] in mask_positions:
            colors.append("red")
        else:
            colors.append("blue")
        frequency_data.append(pos[1])

    x = np.arange(len(position_data))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()

  
    rects1 = ax.bar(x - width/2, frequency_data,  color=colors, width=width, \
        edgecolor='white', label='frequency')


    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Number of Pairs')
    ax.set_xlabel("iSNV Position")
    ax.set_xticks(x)
    ax.set_xticklabels(position_data)
    ax.tick_params(axis = 'x', rotation = 50)

    
    maxn = len(position_data)
    fig.tight_layout() 
    plt.title("Positions")
    plt.xlim(-0.5,maxn-0.5)
    plt.tight_layout()
    #plt.savefig('%s_%s_allele_freq.png'%(donor_filename, recipient_filename), dpi=300, bbox_inches='tight')
    plt.show()
    
#shared_positions('mason_data/', 'default_mask.txt', min_read_depth=10)

masked_shared_variants(shared_variants)
