"""Module containing classes and functions for parsing different types of data"""

#Third party import
import vcf

def _is_valid_lfv(min_read_depth, max_AF, var_reads, total_reads):
    """
    Boolean function called as a helper for determining if a variant fits the
    user-defined parameters for parsing the data.
    """
    # Calculate allele frequency
    freq = var_reads / total_reads

    #If the allele passes restrictions, return True
    if var_reads >= min_read_depth and freq < max_AF:
        return True

    #Otherwise, return False
    return False

class BiallelicData:
    """Dict-like object that stores biallelic data."""

    def __init__(self):
        """Initialize object"""
        self.__dict = dict()

    def __str__(self):
        """Initialize object representation"""
        return str(self.__dict)

    def dict__(self):
        return dict(self)

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

    def store(self, pos, var, freq, var_depth):
        """Stores biallelic data in object"""
        if pos in self.__dict:
            variant = list(self.__dict[pos].values())
            old_freq = variant[0][0]
            if old_freq < freq:
                self.__dict[pos] = {var: [freq, var_depth]}
        else:
            self.__dict[pos] = {var: [freq, var_depth]}

class MultiallelicData:
    """Dict-like object that stores biallelic data."""
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

    def get(self, item):
        """Gets an item stored in object if it is there"""
        return self.__dict.get(item, None)

    def positions(self):
        """Returns a list of all of the positions stored in the object"""
        return self.__dict.keys()

    def store(self, pos, var, freq, var_depth):
        """Stores multiallelic data in object"""
        if pos in self.__dict:
            self.__dict[pos][var] = [freq, var_depth]
        else:
            self.__dict[pos] = {var: [freq, var_depth]}

def extract_lfv(filename, min_read_depth=0, max_AF=1, parse_type="biallelic", store_reference=True, masks=None, mask_status='hide'):
    """
    Extracts variant data from VCF and creates a dictionary storing data
    in the form: {position: {variant: [frequency, depth}}.
    """

    #### Handle Errors #####
    PARSE_TYPES = {"biallelic", "multiallelic"}
    MASK_TYPES = {"hide", "highlight"}
    if min_read_depth < 0 or int(max_AF) > 1 or parse_type not in PARSE_TYPES or str(mask_status) not in MASK_TYPES:
        raise ValueError("Invalid input.")

    #########################

    #Parse mask file if mask file is inputted
    mask_positions = []
    if masks != None:
        with open(masks, "r") as mask_file:
            for line in mask_file:
                comma_split = line.split(",")
                for item in comma_split:
                    nums = item.split("-")
                    if len(nums) == 1:
                        mask_positions.append(int(nums))
                    else:
                        for num in range(nums[0], nums[1]):
                            mask_positions.append(int(num))

    lfv_data = BiallelicData() if parse_type == "biallelic" else MultiallelicData()
    data = vcf.Reader(open(filename, 'r'))
    ref_data = {}

    for record in data:

        # Calculate amount of reads support variant
        var_depth = record.INFO['DP4'][2] + record.INFO['DP4'][3]

        # Calculate amount of reads at particular location
        raw_depth = record.INFO['DP']

        # Find data in VCF file
        pos, var, freq = record.POS, str(record.ALT[0]), float(var_depth / raw_depth)

        #doesn't include masked positions based on user settings
        if pos in mask_positions and mask_status == 'hide':
            break

        # If variant passes restrictions, store data
        if _is_valid_lfv(min_read_depth, max_AF, var_depth, raw_depth):
            lfv_data.store(pos, var, freq, var_depth)

        if store_reference and not pos in ref_data:
            # Calculate how many reads support reference
            ref_depth = record.INFO['DP4'][0] + record.INFO['DP4'][1]

            #Find reference allele
            ref = str(record.REF[0])
            
            # If ref allele passes restrictions, store the data
            if _is_valid_lfv(min_read_depth, max_AF, ref_depth, raw_depth):
                ref_data[pos] = {ref: [(ref_depth / raw_depth), ref_depth]}

   
    # After parsing is complete, make object into a dictionary
    lfv_data = dict(lfv_data)

    

    # If we collected reference data, update lfv_data
    if store_reference:
        for pos in ref_data:
            if pos in lfv_data:
                print("if")
                lfv_data[pos] = ref_data[pos]
            else:
                print("else")
                lfv_data.update(ref_data[pos])

    return lfv_data, mask_positions

def bb_input_data(donor, recip, min_read_depth=0, max_AF=1, parse_type="biallelic", store_reference=True, weighted=False):
    """
    Stores info from parsing VCF files to dictionary.
    """
    donor_data, donor_masks = extract_lfv(donor, min_read_depth=min_read_depth, max_AF=max_AF, parse_type=parse_type, store_reference=store_reference)
    recip_data, recip_masks = extract_lfv(recip, min_read_depth=min_read_depth, max_AF=max_AF, parse_type=parse_type, store_reference=store_reference)

    shared_count = 0

    # Stored as {pos: {var: [donor freq., recip. freq]}} bc Maria had two bb input files 
    # and one required all this info, might change later tho
    bb_data = {} 

    # Iterate through each variant at each position
    for pos in donor_data: 

        bb_data[pos] = {}

        for var in donor_data[pos]:

            # Store donor data
            donor_freq = donor_data[pos][var][0]
            bb_data[pos] = {var: [donor_freq, 0.0]}
            if weighted == True:
                donor_depth = donor_data[pos][var][1]
                bb_data[pos] = {var: [donor_freq, 0.0, donor_depth, 0.0]}

            # If recipient has same variant at same location, store it
            if pos in recip_data and var in recip_data[pos]:
                recip_freq = recip_data[pos][var][0]
                bb_data[pos][var][1] = recip_freq
                shared_count += 1
                if weighted == True:
                    recip_depth = recip_data[pos][var][1]
                    bb_data[pos][var][3] = recip_depth


    return (bb_data, shared_count, donor_masks, recip_masks)
