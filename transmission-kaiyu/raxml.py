import os
import subprocess


class Raxml:
    """
    Wrappers for raxml terminal commands to make running raxml more simple for users.
    Also includes simple methods for getting specific files from raxml outputs.
    """
    def __init__(self, data):
        self.data = data

    def run(self, data):
        """
        Runs raxml on given file and stores files in folder
        """
        #create folder raxml-trees
        os.mkdir("raxml-trees")

        #runs raxml and stores data in raxml-trees
        pass





