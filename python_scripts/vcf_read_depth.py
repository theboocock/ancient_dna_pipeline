#!/usr/bin/env python
#
# Script to take both coverage and fasta files and to mask the reads that are below a certain coverage.`
#
#
#
###
#
# Basic run - take all_sites.vcf, require a minimum read-depth 
#
#
#

import argparse




def main():
    parser = argparser.ArgumentParser(description="Take VCF that contain all site, and output fasta removing all sites with coverage greater than a specified value")
    parser.add_argument("",dest="")
    
if __name__=="__main__":
    main()
