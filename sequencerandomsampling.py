#!/usr/bin/python
# coding: utf-8

#Copyright (c) 2014, Shuji SUZUKI
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#* Redistributions of source code must retain the above copyright notice, this
#  list of conditions and the following disclaimer.
#
#* Redistributions in binary form must reproduce the above copyright notice,
#  this list of conditions and the following disclaimer in the documentation
#  and/or other materials provided with the distribution.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
    parser.add_option("-i", "--input_filenames", dest="input_filenames",
                  help="input filename (if paired files, use \",\" to separate)", metavar="FILE")
    
    parser.add_option("-o", "--output_filename", dest="output_filename",
                  help="output filename", metavar="FILE")
    
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


def WriteSequences(input_handles, file_format, number_sequence, number_sampling_sequences, output_handles):
    sampling_sequence_ids = random.sample(xrange(0, number_sequence), number_sampling_sequences)
    sampling_sequence_ids.sort()

    print sampling_sequence_ids

    for input_handles_i in range(0, len(input_handles)):
        sequences_i = 0
        sampling_sequence_ids_i = 0
        for record in SeqIO.parse(input_handles[input_handles_i], file_format):
            if sequences_i == sampling_sequence_ids[sampling_sequence_ids_i]:
                SeqIO.write(record, output_handles[input_handles_i], file_format)
                sampling_sequence_ids_i = sampling_sequence_ids_i + 1
                if sampling_sequence_ids_i >= number_sampling_sequences:
                    break
            
            sequences_i = sequences_i + 1
    
    return


def main():
    (options, args) = get_options()
    if options.seed == 0:
        random.seed(None)
    else:
        random.seed(options.seed)
    input_filenames = options.input_filenames.split(",")
    number_sampling_sequences = 0
    input_handles = []
    input_handles.append(open(input_filenames[0]))
    number_sequence = GetNumberOfSequences(input_handles[0], options.file_format)
    input_handles[0].close()
    if options.percent == 0.0 and options.number == 0:
        sys.exit("error: unset -p and -n.")
    elif options.percent > 0:
        number_sampling_sequences = int(number_sequence*options.percent/100)
    else:
        number_sampling_sequences = options.number
    
    if number_sampling_sequences > number_sequence:
        sys.exit("error: number sampling sequences is over.")
    input_handles[0] = open(input_filenames[0])
    output_handles = []
    if len(input_filenames) == 2:
        input_handles.append(open(input_filenames[1]))
        output_handles.append(open(options.output_filename + ".1", "w"))
        output_handles.append(open(options.output_filename + ".2", "w"))
    else:
        output_handles.append(open(options.output_filename, "w"))
    
    
    WriteSequences(input_handles, options.file_format, number_sequence, number_sampling_sequences, output_handles)
    
    for i in range(0, len(input_handles)):
        input_handles[i].close()
        output_handles[i].close()


if __name__ == '__main__':
    main()
        
