import os
import subprocess
import glob

class RAxML:
    """
    Wrapper for raxml terminal commands to make running raxml more simple for users.
    Also includes simple methods for getting specific files from raxml outputs.
    """

    def __init__(self, filename, label='raxml', folder='raxml-trees'):
        """
        Initializing attributes for class.
        """
        self.folder = folder
        self.label = label
        self._filename = filename
        self._outputs = {"best": "RAxML_bestTree" + "." + label,\
                        "parsimony": "RAxML_parsimonyTree" + "." + label,\
                        "info": "RAxML_info" + "." + label}
        
    def getPath(self, item):
        """
        Returns path specified tree file using keyword "item".
        """
        path = os.path.join(os.getcwd(), self.folder, self._outputs[item])
        if os.path.exists(path):
            return path
        raise ValueError("Cannot find file. Try using run() method to create RAxML files first.")

    def run(self, overwrite=False): 
        """
        Runs raxml on given file and stores files in folder.
        """
        path1 = os.path.join(os.getcwd(), self.folder)
        if overwrite:
            if os.path.exists(path1):
                filtered_dir = [f for f in os.listdir(path1) if f.endswith("."+self.label)]
                for files in filtered_dir:
                        os.remove(os.path.join(path1, files))
            else:
                raise ValueError("Directory is missing. Try running RAxML first with run() method.")
        else:
            os.mkdir(self.folder)

        #locate absolute path of file (raxml only takes absolute paths)
        path2 = os.path.join(os.getcwd(), self._filename)
        
        #runs raxml and stores data in raxml-trees
        command = f"raxmlHPC -s {path2} -w {os.path.join(os.getcwd(), self.folder)}\
                 -n {self.label} -m GTRGAMMA -p 20"
        subprocess.run(command.split())
        
        #print out command
        msg = f'Ran RAxML on {self._filename} and stored files in directory: {self.folder}' + '.'
        print(msg)
        


