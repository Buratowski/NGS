#python 3.4
#Dependencies:
#   numpy
#   matplotlib
#   h5py
#
#Collection of scripts to load and analyze NGS data

import numpy as np
import matplotlib.pyplot as plt
import h5py

#Create list of tuples containing chromosome names and sizes (chromosomes)
#Values are for S.cerevisiae

list_of_chromosomes = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV',
                          'XVI', 'MT'] # The names of all chromosomes
size_of_chromosomes = [230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816,
                           1078177, 924431, 784333, 1091291, 948066, 85779] #the sizes of all chromosomes
chromosomes=[(chrom,size) for chrom,size in zip(list_of_chromosomes,size_of_chromosomes)]

def load_track(BED):
    """ This function creates a dictionary of arrays and populates it with data from the BED file
    Arguments:
          BED string of BDG file name (file should not contain header)
    Returns:
          Dictionary of numpy arrays with chromosome sizes and position values based on BDG file data
    Example of usage:
          track_data=load_track(<BDG file>)
    """
    
    data=dict()
    for item in chromosomes:
        data[item[0]]=np.zeros(item[1]) #create an empty array with the size of the chromosome
    with open(BED) as data_file:
        for line in data_file:
            row=line.split('\t')
            data[row[0]][int(row[1]):int(row[2])]=float(row[3]) #populate array
    return data
