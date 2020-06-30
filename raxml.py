import os
import subprocess
import shutil

"""
A module containing functions relating to RAxML.
"""

def run(filename, label='raxml', outputDir='raxml-trees', customCMD = ''): 
    """
    A function for running RAxMl. Writes outputs to outputDir.

    Args:
        filename (str): the name of the file that RAxMl will analyze
        label (str, optional): The label tagged at the end of the output files. Defaults to 'raxml'.
        outputDir (str, optional): Name of the folder in which RAxML will output to. Defaults to 'raxml-trees'.
        customCMD (str, optional): A custom RAxML command. Defaults to ''.
    """

    # Remove directory if it exists
    if os.path.exists(outputDir):
        shutil.rmtree(outputDir)
    os.mkdir(outputDir)

    if customCMD:
        # If user specifies command, then we run that in command line
        cmnd = customCMD
    
    else:
        # Otherwise, run RaxML with default arguments
        path = os.path.join(os.getcwd(), filename)
        cmnd = f"raxmlHPC -s {path} -w {os.path.join(os.getcwd(), outputDir)} -n {label} -m GTRGAMMA -p 20"

    subprocess.run(cmnd.split())
    
    # Print out when finished
    msg = f'Ran RAxML on {filename} and stored files in directory: {outputDir}' + '.'
    print(msg)
        
