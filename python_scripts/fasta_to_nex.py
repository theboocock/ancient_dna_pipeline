#!/usr/bin/env python
#
##
#
# Remove interleaved state from a fasta sequences
#
# @author James Boocock



import argparse


def fasta_to_nexus(

def main():
    parser = argparse.ArgumentParser(description='Un interleave_nexus')
    
    parser.add_argument('-i','--input',dest='input_fasta',help='input fastafile')
    parser.add_argument('-o','--output',dest='output_nex',help='output nexus file')
    args = parser.parse_args()
    load_sequences(args.input_fasta)
    fasta_to_nexus(args.input_fasta, args.output)

