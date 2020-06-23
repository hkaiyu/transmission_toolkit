
def shared_variants(lfv_list1, lfv_list2):
    """
    Inputs:
        lfv_dict1: a list of tuples containing low frequency variants in the form (position, variant, frequency)
        lfv_dict2: a list of tuples containing low frequency variants in the form (position, variant, frequency)
    Outputs: 
        shared: a list of tuples with all of the shared variants between the two input files in the form (position, variant, frequency1, frequency2)
    """
    #adds filenames as first two elements of shared
    shared = [lfv_list1[0], lfv_list2[0]]
    #get list of shared variants
    for num1 in range(1, len(lfv_list1) - 1):
        for num2 in range(1, len(lfv_list2) - 1):
            #add to list of shared variants
            if lfv_list1[num1][0] == lfv_list2[num2][0] and lfv_list1[num1][1] == lfv_list2[num2][1]:
                shared.append((lfv_list1[num1][0], lfv_list1[num1][1], lfv_list1[num1][2], lfv_list2[num2][2]))
    return shared