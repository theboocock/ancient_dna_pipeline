#!/usr/bin/env python
#
# Plot exponential decay data read bam file. 
#

import pysam
import argparse

def main():
    args = argparse.ArgumentParser(description="Exponential decay Script for a bam file")
    args.add_argument('bam',help="Bam input file")
    exponential_decay

if __name__=="__main__":
    main():
