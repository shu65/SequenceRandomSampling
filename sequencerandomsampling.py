#!/usr/bin/python
# coding: utf-8

import sys
import random
from optparse import OptionParser

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
from numpy.core.fromnumeric import sort


def get_options():
    parser = OptionParser()
    parser.add_option("-i", "--input_filename", dest="input_filename",
                  help="input fasta filename", metavar="FILE")
    
    parser.add_option("-o", "--output_filename", dest="output_filename",
                  help="Output fasta filename", metavar="FILE")
    
    parser.add_option("-p", "--percent", dest="percent",
                  help="sampling percent (%)", type="float", default=0.0)
    
    parser.add_option("-n", "--number", dest="number",
                  help="number of subsequences ", type="int", default=0)
    
    parser.add_option("-f", "--sequence_file_format", dest="file_format",
                  help="sequence file format (fasta, fastq)", default="fasta")
    
    parser.add_option("-s", "--seed", dest="seed",
                       help="random seed", type="int", default=0)
    
    return parser.parse_args()


def GetNumberOfSequences(input_handle, file_format):
    return len([None for record in SeqIO.parse(input_handle, file_format)])


def WriteSequences(input_handle, file_format, number_sequence, number_sampling_sequences, output_handle):
    sampling_sequence_ids = random.sample(xrange(0, number_sequence), number_sampling_sequences)
    sampling_sequence_ids.sort()

    print sampling_sequence_ids

    sequences_i = 0
    sampling_sequence_ids_i = 0
    for record in SeqIO.parse(input_handle, file_format):
        if sequences_i == sampling_sequence_ids[sampling_sequence_ids_i]:
            SeqIO.write(record, output_handle, file_format)
            sampling_sequence_ids_i = sampling_sequence_ids_i + 1
            if sampling_sequence_ids_i >= number_sampling_sequences:
                return
            
        sequences_i = sequences_i + 1
    
    return


def main():
    (options, args) = get_options()
    if options.seed == 0:
        random.seed(None)
    else:
        random.seed(options.seed)
    number_sampling_sequences = 0
    input_handle = open(options.input_filename)
    number_sequence = GetNumberOfSequences(input_handle, options.file_format)
    input_handle.close()
    if options.percent == 0.0 and options.number == 0:
        sys.exit("error: unset -p and -n.")
    elif options.percent > 0:
        number_sampling_sequences = int(number_sequence*options.percent/100)
    else:
        number_sampling_sequences = options.number
    
    if number_sampling_sequences > number_sequence:
        sys.exit("error: number sampling sequences is over.")
    input_handle = open(options.input_filename)
    output_handle = open(options.output_filename, "w")
    
    WriteSequences(input_handle, options.file_format, number_sequence, number_sampling_sequences, output_handle)
    
    input_handle.close()
    output_handle.close()


if __name__ == '__main__':
    main()
        
