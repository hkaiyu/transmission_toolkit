import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def bar_plot_input(donor, recipient):
    """
    Inputs:
        donor: a list of tuples containing low frequency variants in the form (position, variant, frequency) with the first element being the file name
        recipient: a list of tuples containing low frequency variants in the form (position, variant, frequency) with the first element being the file name
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

def bar_plot(bar_plot_input):
    """
    Inputs:
    a list of tuples of donor variants and their frequencies in the donor and recipient in the form (position, variant, donor frequency, recipient frequency2)
    Outputs:
    a bar plot containing the donor variants and frequencies in the donor and recipient 
    """
    
    filename1 = bar_plot_input[0]
    filename2 = bar_plot_input[1]

    positions = []
    frequency1 = []
    frequency2 = []
    allele = []

    for num in range(2, len(bar_plot_input)):
        positions.append((bar_plot_input[num][0], bar_plot_input[num][1]))
        allele.append(bar_plot_input[num][1])
        frequency1.append(round(bar_plot_input[num][2], 2))
        frequency2.append(round(bar_plot_input[num][3], 2))

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