"""Utility functions and classes"""

import ntpath

def get_seq(sequence):
    """
    Returns string representation of sequence genome given a FASTA file.
    Assumes only one sequence in the file
    """    
    consensus = str()

    with open(sequence, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
        for line in lines[1:]:
            consensus += line
            if line and line[0] == '>':
                raise ValueError('File should only include one sequence.')

    return consensus
    
def getpathleaf(path):
    '''
    Returns the leaf of a given path.
    For example, if inputted, home/user/.../file, it
    will return file
    '''
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)