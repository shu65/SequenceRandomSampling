#!/usr/bin/python
# coding: utf-8

import sys
from optparse import OptionParser

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC

def get_options():
    parser = OptionParser()
    parser.add_option("-a", "--ssearch_result", dest="ssearch_result_filename",
                  help="SSEARCH result file name", metavar="FILE")
    
    parser.add_option("-c", "--comparing_result", dest="comparing_result_filename",
                  help="Result file name to compare", metavar="FILE")
    
    parser.add_option("-t", "--tool", dest="tool",
                  help="comparing tool name : blast, ghostx", type="string")
    
    parser.add_option("-x", "--x_axis", dest="x_axis",
                  help="x axis : E-value (e), Bit-score (b), Identity (i)", type="string")
    
    parser.add_option("-s", "--x_start", dest="x_start", default= -1.0,
                  help="x start value", type="float")
    
    parser.add_option("-e", "--x_end", dest="x_end", default= -1.0,
                  help="x end value", type="float")
    
    parser.add_option("-p", "--print_mishit", action="store_true", dest="print_mishit",
                  help="print mis hit result")

    return parser.parse_args()

def main():
    (options, args) = get_options()
    ssearch_file = open(options.ssearch_result_filename)
    ssearch_reader = SSearchResultReader(ssearch_file)
    x_axis = ""
    if options.x_axis == "e":
        x_axis = "e_value"
    elif options.x_axis == "b":
        x_axis = "bit_score"
    elif options.x_axis == "i":
        x_axis = "identity"
    else:
        sys.exit("errors : x axis is invalided")
    x_start, x_end, iterating_op, iterating_value = get_x_axis_default(ssearch_reader, x_axis)
    ssearch_file.close()
    
    if options.x_start != -1:
        x_start = options.x_start
    if options.x_end != -1:
        x_end = options.x_end
    
    if x_start > x_end:
        sys.exit("errors : x start and end")
        
    x_range = get_x_range(x_axis, iterating_op, iterating_value, x_start, x_end)
    
    ssearch_file = open(options.ssearch_result_filename)
    ssearch_reader = SSearchResultReader(ssearch_file)
    comparing_file = open(options.comparing_result_filename)
    comparing_reader = BlastResultReader(comparing_file);
    # comparing_reader = ExtendingBlastResultReader(comparing_file);
    hit_list, sum_list, number_ignored_alignments = compare_ssearch(ssearch_reader, comparing_reader, x_axis, x_range, options.print_mishit)
    output(x_range, hit_list, sum_list, number_ignored_alignments)
    
    ssearch_file.close()
    comparing_file.close()
    

usage = "usage: %prog [options]"
parser = OptionParser(usage)

parser.add_option(
    "-i", "--input", 
    action="store",
    type="string",
    dest="input",
    help="input fasta file"
)


parser.add_option(
    "-o", "--output",
    action="store",
    type="string", 
    dest="output",
    help="output fasta file prifix",
)


parser.add_option(
    "-n", 
    action="store",
    type="int",
    dest="number",
    help="division number"
)

(options, args) = parser.parse_args()

input_file  = options.input
output_file = options.output + "_"

number_records = 0
input_handle = open(input_file, "r")
for record in SeqIO.parse(input_handle, "fasta"):
    number_records = number_records + 1

stride = int((number_records + options.number - 1)/options.number)

input_handle.close()
input_handle = open(input_file, "r")
for i in range(options.number):
    records = []
    chunk_size = 0
    for record in SeqIO.parse(input_handle, "fasta"):
        records.append(record)
        chunk_size = chunk_size + 1
        if chunk_size >= stride :
            break

    out_handle = open(output_file + str(i), "w")
    for record in records[:] :
        SeqIO.write(record, out_handle, "fasta")

    out_handle.close()

input_handle.close()
'''
records = []
for record in SeqIO.parse(input_file, "fasta"):
    records.append(record)

number_records = len(records)
stride = int((number_records + options.number - 1)/options.number)

for i in range(options.number):
    out_handle = open(output_file + str(i), "w")
    for record in records[stride*i:stride*(i + 1)] :
        SeqIO.write(record, out_handle, "fasta")

    out_handle.close()


Created on 2013/10/07

@author: shu
'''

        
if __name__ == "__main__":
    main()
        
        