"""Utility functions and classes"""

class Data:
    """Class for storing data in dict-like object"""
    def __init__(self):
        """Initialize object"""
        self.__dict = dict()

    def __str__(self):
        """Initialize object representation"""
        return str(self.__dict)

    def __getitem__(self, key):
        """Returns item at key in object"""
        return self.__dict[key]

    def __setitem__(self, key, value):
        """Sets key-value pair in object"""
        self.__dict[key] = value

    def __iter__(self):
        """Returns iterable of object"""
        return iter(self.__dict.items())

    def __len__(self):
        """Returns length of object"""
        return len(self.__dict)

    def __delitem__(self, key):
        """Deletes item from object"""
        del self.__dict[key]

    def positions(self):
        """Returns a list of all of the positions stored in the object"""
        return self.__dict.keys()

class Biallelic(Data):
    """Dict-like object that stores biallelic data."""

    def __init__(self):
        """Initialize object"""
        super().__init__()
        self.__dict = dict()

    def store(self, pos, var, freq, var_depth):
        """Stores biallelic data in object"""
        if pos in self.__dict:
            variant = list(self.__dict[pos].values())
            old_freq = variant[0][0]
            if old_freq < freq:
                self.__dict[pos] = {var: [freq, var_depth]}
        else:
            self.__dict[pos] = {var: [freq, var_depth]}

class Multiallelic(Data):
    """Dict-like object that stores biallelic data."""
    def __init__(self):
        """Initialize object"""
        super().__init__()
        self.__dict = dict()

    def store(self, pos, var, freq, var_depth):
        """Stores multiallelic data in object"""
        if pos in self.__dict:
            self.__dict[pos][var] = [freq, var_depth]
        else:
            self.__dict[pos] = {var: [freq, var_depth]}
