import matplotlib.pyplot as plt
import datetime
import numpy as np
from parsers import parse_biallelic as biallelic, parse_multiallelic as multiallelic

def make_bar_plot(donor_filename, recipient_filename, parse_type, min_AF = 0.03, max_AF = 0.5):
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
    vcf_data = parse_type(donor_filename, recipient_filename)[0]
    positions, donor_freqs, recipient_freqs = [], [], []

    for pos in sorted(vcf_data):
        for nuc in vcf_data[pos]:
            positions.append((pos, nuc))
            donor_freqs.append(round(vcf_data[pos][nuc][0], 3))
            recipient_freqs.append(round(vcf_data[pos][nuc][1], 3))

    x = np.arange(len(positions))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, donor_freqs,  color='#7f6d5f', width=width, \
        edgecolor='white', label=donor_filename)
    rects2 = ax.bar(x + width/2, recipient_freqs, color='#557f2d', width=width, \
        edgecolor='white', label=recipient_filename)
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Allele Frequencies')
    ax.set_title('Allele Frequency vs. Position')
    ax.set_xticks(x)
    ax.set_xticklabels(positions)
    ax.tick_params(axis = 'x', rotation = 35)
    ax.legend()


    def autolabel(rects):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)

    fig.tight_layout()
    
    plt.savefig('%s_%s_allele_freq.png'%(donor_filename, recipient_filename), dpi=fig.dpi, bbox_inches='tight')
    plt.show()
