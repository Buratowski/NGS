#python3.3
#Executable script to remove 5' bases from all reads in fastq file
#Running instructions
#   fastq_trimmer -f <number of first base to include, usually lenght of barcode +2> -i <input file, fastq> -o <output file, fastq>

import argparse

parser = argparse.ArgumentParser(description='Split FASTQ file according to barcodes')
parser.add_argument('-f', dest='first', help='first base to keep',nargs=1,required=True)
parser.add_argument('-i', dest='input_file', help='input file',nargs=1,required=True)
parser.add_argument('-o', dest='output_file', help='output file prefix',nargs=1,required=True)
args = parser.parse_args()

# Parse arguments from command line execution
first=int(args.first[0])
input_file=args.input_file[0]
output_file=args.output_file[0]


def create_output_files(input_file,output_file,first):
    """
    Creates new 5' trimmed fastq file, removes both bases and corresponding qualities:
    Arguments:
              input_file: string of input file name
              output_file: string to rename output file, will overwrite if exists
              first: int of first base to include in reads in output file
    Returns:
              None
    """
    
    starting_file=open(input_file,'r')
    out_file=open(output_file,'w')
    header='start'
    while header:
        header=starting_file.readline()
        sequence=starting_file.readline()
        header=header+sequence[first:]
        header=header+starting_file.readline()
        header=header+starting_file.readline()[first:]
        out_file.write(header)
    out_file.close()

if __name__=='__main__':
    create_output_files(input_file,output_file,first)
