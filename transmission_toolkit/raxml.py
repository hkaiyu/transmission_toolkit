"""
A module containing functions for running RAxML
"""
#Standard library imports
import os
import subprocess
import glob

def run_raxml(filename, label='raxml', output_dir='raxml-trees', custom_cmd=''):
    """
    A function for running RAxMl. Writes outputs to output_dir.

    Args:
        filename (str): the name of the file that RAxMl will analyze
        label (str, optional): The label tagged at the end of the output files.
        output_dir (str, optional): Name of the folder in which RAxML will output to.
        custom_cmd (str, optional): A custom RAxML command.
    """

    # If directory already exists, delete all files with same label
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
        path = os.path.join(os.getcwd(), filename)
        cmnd = f"raxmlHPC -s {path} -w {os.path.join(os.getcwd(), output_dir)} \
            -n {label} -m GTRGAMMA -p 20"
    subprocess.run(cmnd.split(), check=True)

    # Print out when finished
    msg = f'Ran RAxML on {filename} and stored files in directory: {output_dir}' + '.'
    print(msg)
