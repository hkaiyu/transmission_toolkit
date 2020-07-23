"""Module for extracting data from fasta files and aligning seqences"""
import os
import glob
import errno
import subprocess
from pathlib import Path

from transmission_toolkit.utils import getpathleaf
from transmission_toolkit.VCFtools import extract_lfv, build_consensus

FATSA_EXTENSIONS = {'.fasta', '.fna', '.ffn', '.faa', '.frn'}
LINE_WRAP = 80 # Max. length of each line in fasta file 
THREADS = 2

def write_fasta(vcf_path, reference, output_dir="", line_length=LINE_WRAP, **filters):
    """
    Writes a FASTA file with a VCF file and reference FASTA file as input.
    """
    # If user specifies output directory
    if output_dir:

        # Check if directory already exists
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

    name = getpathleaf(vcf_path).split('.')[0]
    path = os.path.join(output_dir, name + '.fna')
    seq = build_consensus(vcf_path, reference, **filters)

    with open(path, 'w') as f:
        f.write(f'>{name}\n')
        for i in range(0, len(seq), line_length):
            f.write(seq[i: i+line_length] + '\n')

'''
path = os.path.join('example_data/mason_data')
for f in os.listdir(path):
    fn = os.path.join(path, f)
    write_fasta(fn, 'example_data/sequence.fasta', output_dir='example_data/tmpdata')
'''

class FastaAligner:
    """
    Class for aligning FASTA files.

    This class is meant for handling FASTA files with one sequence per file, and will raise
    an error if given a multi-FASTA file.
    """
    def __init__(self, fasta_dir):
        if os.path.isdir(fasta_dir):
            self.dir = fasta_dir
        else:
            raise ValueError("The specified path should be a directory of fasta files.")

    def align(self, reference, output_dir='', threads=THREADS):
        """
        Given a reference genome file, aligns all fasta files in directory using Parsnp.
        """
        if not os.path.exists(reference):
            raise FileNotFoundError(f"{reference}")
        if output_dir and not os.path.exists(output_dir):
            os.mkdir(output_dir)
        
        output = os.path.join(os.getcwd(), output_dir)
        cmnd = f"parsnp -d {self.dir} -r {reference} -o {output} -p {threads}"
        subprocess.call(cmnd.split())

        path2xmfa = os.path.join(output, 'parsnp.xmfa')
        cmnd2 = f"harvesttools -x {path2xmfa} -M " + os.path.join(output_dir,'parsnp.mfa')
        subprocess.call(cmnd2.split())

class FastaRecord:
    """
    Simple class for storing and accessing data from a FASTA record.
    """
    def __init__(self, name, sequence):
        if not isinstance(name, str):
            raise TypeError('Name should be a string.')
        if not isinstance(sequence, str):
            raise TypeError('Sequence should be a string.')

        self.name = name
        self.seq = sequence
    
    def set_id(self, identifier):
        if not isinstance(identifier, (str, int, float)):
            raise TypeError('ID should be a string, integer, or float.')
        self.id = identifier

class MultiFastaParser:
    """
    Makes it easier to access data in a multi-FASTA file.
    """
    def __init__(self, multifasta):
        if os.path.isfile(multifasta):
            self.fasta = multifasta
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), multifasta)
        
        self.records = []

        # Populates self.records with FastaRecord objects
        with open(self.fasta, 'r') as f:
            lines = [line.strip() for line in f.readlines()]
            length = len(lines)
            for i in range(0, length - 1):
                if lines[i][0] == '>':
                    j = 1
                    while i+j < length:
                        if lines[i+j][0] == '>':
                            break
                        j += 1
                    consensus = ''.join(line for line in lines[i+1:i+j])
                    name = lines[i][1:]
                    self.records.append(FastaRecord(name, consensus))

    def get_groups(self):
        """
        Method that parses multi-FASTA file and groups records by sequence in dictionary.
        """
        groups = dict()
        for record in self.records:
            if record.seq in groups:
                groups[record.seq].add(record.name.split('.')[0])
            else:
                groups[record.seq] = {record.name.split('.')[0]}
        return groups.values()

    def infer_phylogeny(self, output_dir='', label='tree', threads=THREADS, custom_cmd=''):
        """
        Runs RAxML on the multi-fasta file and stores output in output_dir.
        """
        # If directory already exists, delete all files with same label
        if output_dir:
            if os.path.exists(output_dir):
                for fname in glob.glob(os.path.join(output_dir, f"*.{label}")):
                    os.remove(os.path.join(output_dir, fname))
            else:
                os.mkdir(output_dir)

        # If user specifies command, then we run that in command line
        if custom_cmd:
            cmnd = custom_cmd

        # Otherwise, run RaxML with default arguments
        else:
            cwd = os.getcwd()
            path = os.path.join(cwd, self.fasta)
            cmnd = f"raxmlHPC -s {path} -w {os.path.join(cwd, output_dir)} \
                -n {label} -m GTRGAMMA -p {threads}"
        subprocess.call(cmnd.split())

        # Print out when finished
        msg = f'Ran RAxML on {self.fasta} and stored files in directory: {output_dir}' + '.'
        print(msg)

#data = MultiFastaParser('TransmissionViz/Parsnp/parsnp.mfa')
