#!/usr/bin/env python
#
# Plot exponential decay data read bam file. 
#

import pysam
import argparse

def main():
    parser = argparse.ArgumentParser(description="Exponential decay Script for a bam file")
    parser.add_argument('bam',help="Bam input file")
    args = parser.parse_args()
    
    exponential_decay()

if __name__=="__main__":
    main():
