import numpy as np
import argparse
import time
import math
import sys
from Bio import SeqIO
#import random
#import sa

def read_seqfile(input_file):
    sys.stderr.write("\tReading Sequence File {0} ... ".format(input_file))
  
    records = list(SeqIO.parse(input_file, "fasta"));
    nseq = len(records);
    sys.stderr.write(" Number of sequences: %d\n"%nseq)

    d_arr = np.empty(nseq, dtype='object')
    seq_labels = np.empty(nseq, dtype='object')
    i = 0;
    for record in SeqIO.parse(input_file, "fasta"):
        d_arr[i] = record.seq;
        seq_labels[i] = record.id
        i += 1;
    return (d_arr,seq_labels)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seqfile', help='Inupt filename')
    parser.add_argument('bowtiefile', help='Bowtie seed map file')
    #parser.add_argument('nSeq', help='Number of sequences', default=100, type=int)
    #parser.add_argument('nOverlap', help='Length of overlap', type=int)
    parser.add_argument('outpath', help='Output file to store overlaps')
    args = parser.parse_args()

    global seqarray
    global seqlabels
    seqarray,seqlabels = read_seqfile(args.seqfile);
    
    seqlbl_dictionary = {};
    for i in range(seqarray.shape[0]):
        seqlbl_dictionary[seqlabels[i]] =i;
    
    seed_file = open(args.outpath+"/seed_reads.fa","w")
    with open(args.bowtiefile) as f:
        for line in f:
            #print(line);
            line_content = line.split("\t");
            #print(line_content[0:4]);
            if line_content[0] in seqlbl_dictionary:
                #print('>',line_content[0]);
                seed_file.write(">%s\n"%(line_content[0]))
                #seed_file.write(line_content[0]);
                #seed_file.write("\n");
                seq_index = seqlbl_dictionary[line_content[0]];
                #print(seqarray[seq_index]);
                seed_file.write("%s\n"%(seqarray[seq_index]));
    
    seed_file.close();     


if __name__ == '__main__':
    main()

