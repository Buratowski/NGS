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
          
    Legacy function not needed if using hdf5 track loading
    """
    
    data=dict()
    for item in chromosomes:
        data[item[0]]=np.zeros(item[1]) #create an empty array with the size of the chromosome
    with open(BED) as data_file:
        for line in data_file:
            row=line.split('\t')
            data[row[0]][int(row[1]):int(row[2])]=float(row[3]) #populate array
    return data

hdf5_file=h5py.File("datasets.hdf5","w") #create a hdf5 file

def create_hdf5_track(dataset_name,hdf5_file_handler):
    """Creates hdf5 track containing subgroups of zero valued arrays, one for each chromosome
    requires a previous created dictionary of chromosome names/sizes named chromosomes
    Arguments:
      dataset_name string with the name of dataset to create
      hdf5_file_handler variable pointing to open hdf5 file
     Returns:
      None (modifies opened hdf5 file)
     Usage example:
      create_hdf5_track(track_name,hdf5_file_handler)
    """
    temp=hdf5_file_handler.create_group(dataset_name) #create a first group with the name of our dataset
    for item in chromosomes:
        temp.create_dataset(item[0],(item[1],),dtype='f')
        
def create_hdf5_track_compressed(dataset_name,hdf5_file_handler):
    """Similar to create_hdf5_track only file will be compressed, usefull for large/sparce data
    requires a previous created dictionary of chromosome names/sizes named chromosomes
    Arguments:
      dataset_name string with the name of dataset to create
      hdf5_file_handler variable pointing to open hdf5 file
     Returns:
      None (modifies opened hdf5 file)
     Usage example:
      create_hdf5_track_compressed(track_name,hdf5_file_handler)
    """
    temp=hdf5_file_handler.create_group(dataset_name) #create a first group with the name of our dataset
    for item in chromosomes:
        temp.create_dataset(item[0],(item[1],),dtype='f',compression='gzip', compression_opts=9)
        
def load_track_in_hdf5(BED,track_name,hdf5):
    """
    Function to populate hdf5 file with data from BDG file
    Requires a previous created dictionary of chromosome names/sizes named chromosomes
    Arguments:
      BED string of BDG filename (BDG should not have header)
      track_name string of dataset name in open hdf5 file
      hdf5 hdf5_file_handler variable pointing to open hdf5 file
    Returns:
      None (modifies opened hdf5 file)
    Usage example:
      load_track_in_hdf5(BDF filename,track name,hdf5_file_handler)
    """
    
    data=dict() 
    for item in chromosomes:
        data[item[0]]=np.zeros(item[1]) #create an empty array with the size of the chromosome
    with open(BED) as data_file:
        for line in data_file:
            row=line.split('\t')
            data[row[0]][int(row[1]):int(row[2])]=float(row[3]) #populate array
    for key,value in data.items():
            hdf5[track_name][key][...] = value #populate hdf5 file
        
#Gene class
#Blueprint to create Gene Class instances
class Gene():
    '''The simplest description of a gene should contain
    its name, chromossome, star location, end location and strand'''
    def __init__(self,name,chromosome,start,end,strand):
        self.name=name
        self.chromosome=chromosome
        self.start=int(start)
        self.end=int(end)
        self.strand=strand
    def __repr__(self):
        return "Gene name:{} ,Start:{}, End:{}, Chromosome:{}, Strand{}".format(self.name,self.start,self.end,
                                                                                self.chromosome,self.strand)

def create_gene_list(bed):
    """
    Creates list of GENE instances based on BED file containing names and coordinates
    Arguments:
      bed string of tab delimited file containing fields in this order:
        Chromosome_name start end strand  name
      file should contain header
    Returns:
      list of gene class instances
    """
    gene_list=[] # this will be the gene list we will populate
    with open(bed) as data:
        data.readline() #get rid of the first line since is an header
        for line in data:
            temp=line.split('\t')
            """Two things to notice, first Start is always smaller than End meaning that the names
            dont reflect the biological meaning, this is important because they should be
            interpreted together with the strand value, second the chromosome location uses 
            the chrX notation which is unlike the description in the dataset we are going
            to use which only uses X for chromosome names"""
            gene_list.append(Gene(temp[4].rstrip(),temp[0][3:],temp[1],temp[2],temp[3]))
    return gene_list

def create_matrix(nucleotide_range,genes,data,track):
    """Creates a numpy 2D matrix containg values of data for each gene that can be used to generate heatmaps/anchor plots
      x=nucleotide range
      y=gene
    Arguments:
      nucleotide_range int to retrieve nucleotide range from -nucleotide range to +nucleotide range, 0 is gene start
      genes list of gene class instances
      data hdf5 file handler
      track string name of track
    Returns:
      2D numpy matrix
    Considerations:
      genes for which range falls over the chromosome genes will not be used
      order of row matrix creation relies on gene list order
      faster implementation should create first a 2D matrix to start
    """
    #first we create a simple array with the nucleotide range
    matrix=np.zeros(2*nucleotide_range) #the range extends on both sides
    n=0 #we create a variable to hold the number of genes we are not plotting (see later)
    for item in genes:
        try: #it is here because some range will fall out of the chromosome range
            if item.strand=='+':
                temp_array=np.array(data[track][item.chromosome][item.start-nucleotide_range:item.start+nucleotide_range])
            else:
                temp_array=np.array(data[track][item.chromosome][item.end-nucleotide_range:item.end+nucleotide_range])[::-1]
            #this is probably not the most efficient way to do it but we don't how many genes can be loaded
            matrix=np.vstack((matrix,temp_array))
        except ValueError:
            n+=1
    print("Unable to load {} gene(s).".format(n))
    return matrix[1::]

#To plot heat map with create_matrix return matrix
#Remember that imshow will reverse y order
fig,ax=plt.subplots(figsize=(3,5))
image_plot=ax.imshow(matrix,aspect='auto',interpolation='none')
#cmap can be change to other matplotlib cmaps check http://matplotlib.org/examples/color/colormaps_reference.html
image_plot.set_cmap('Reds')

#To plot anchor plot with create_matrix return matrix
plt.plot(range(0,matrix.shape[1]),matrix.mean(axis=0))
         
      
