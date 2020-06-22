import subprocess
import os
import argparse
import csv
from collections import defaultdict
import vcf
import pysam
import numpy as np
from Bio import SeqIO

def simulate_consensus(sample_prefix, reference_file, output_dir, filtered_vcf_dir, consensus_dir):
    '''simulate consensus based on vcf file'''
    print(sample_prefix)
    # vcf file zipping, indexing, and consensus generation
    subprocess.run([
            "bgzip",
            os.path.join(output_dir, filtered_vcf_dir, f"{sample_prefix}.vcf")],
        check=True)
    subprocess.run([
            "bcftools",
            "index",
            os.path.join(output_dir, filtered_vcf_dir, f"{sample_prefix}.vcf.gz")],
        check=True)
    subprocess.run([
            "bcftools",
            "consensus",
            os.path.join(output_dir, filtered_vcf_dir, sample_prefix + ".vcf.gz"),
            "-f",
            reference_file,
            "-o",
            os.path.join(output_dir, consensus_dir, sample_prefix + ".fasta")],
        stdout=open(os.path.join(output_dir, "simulation.log"), "a"),
        stderr=open(os.path.join(output_dir, "simulation.err"), "a"),
        check=True)
    # unzip vcf file
    subprocess.run([
            "bgzip",
            "-d",
            os.path.join(output_dir, filtered_vcf_dir, f"{sample_prefix}.vcf.gz")],
        check=True)

def variant_filtering(sample_metadata, output_dir, vcf_dir, filtered_dir, min_af_threshold):
    vcf_files = vcf_dir
    filtered_vcf_files = os.path.join(output_dir, filtered_dir)
    if not os.path.exists(filtered_vcf_files):
        os.mkdir(filtered_vcf_files)

    for sample in sample_metadata:
        vcf_input = os.path.join(vcf_files, f"{sample}.vcf")
        if os.path.exists(os.path.join(filtered_vcf_files, f"{sample}.vcf")):
            os.remove(os.path.join(filtered_vcf_files, f"{sample}.vcf"))
        if os.path.exists(vcf_input):
            subprocess.run([
                    "lofreq",
                    "filter",
                    "-a",
                    f"{min_af_threshold}",
                    "-i",
                    vcf_input,
                    "-o",
                    os.path.join(filtered_vcf_files, f"{sample}.vcf")])

def parse_dir(vcf_dir):
    sample_files = os.listdir(vcf_dir)
    samples = []
    for sample_file in sample_files:
        samples.append(sample_file.rstrip(".vcf"))
    return samples

