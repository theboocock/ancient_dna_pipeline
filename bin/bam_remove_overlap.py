#!/usr/bin/env/python
#
#
#  This program performs 2 key tasks for ancient DNA sequencing.
#  Lowers the quality score of reads as in the neanderthal paper.
#  Also returns one file 
#


import argparse





def main():
    parser = argparse.ArgumentParser(description="Filter paired-end reads from fastq for ancient DNA purposes")
    parser.add_argument('-i','--input',dest='vcf_input',description='VCF input file')
    parser.add_argument('-o','--output',dest='vcf_output',description='VCF output')
    parser.add_argument('
