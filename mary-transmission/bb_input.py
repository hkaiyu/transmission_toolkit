import numpy as np

def bottleneck_input(bb_input):
    """
    Inputs:
        bb_input: a list of tuples with all of the donor variants and their frequencies in the recipient in the form 
        (position, variant, donor frequency, recipient frequency), assumes filename as first element
    Output:
        bottleneck.txt: a txt file in the form of a bb_bottleneck input file
    """
    bb = open("bottleneck.txt", "w+")
    for num in range(2, len(bb_input) - 1):
        bb.write(str(bb_input[num][2]) + "\t" + str(bb_input[num][3]) + "\n")
    bb.close()