#!/usr/bin/env python
#
#
# Calculate the percentage of reads mapping to all mtDNA, other than the first one..
#
# Currently works for mtDNA only.

import argparse
import pysam

def get_reads(bam_file, reference_name,mapq=20):
    bamfile = pysam.AlignmentFile(bam_file)
    reads = bamfile.fetch(reference_name)
    count_total_ref = 0
    for read in reads:
        if (read.mapping_quality >= mapq):
            count_total_ref += 1
    head = bamfile.header
    for item in head['SQ']:
        if(item['SN'] != reference_name):
            print item['SN'],
            count_contam = 0
            for read in bamfile.fetch(item['SN']):
                if read.mapping_quality >= mapq:
                    count_contam +=1
            print(count_contam/float(count_total_ref))

def main():
    parser=argparse.ArgumentParser(description="Compares against other reference Genomes for evidence of Contamination")
    parser.add_argument("bam_file",help="Bam file input, contains the multiple references")
    parser.add_argument("-r","--reference",dest="reference_name",help="Reference_name")
    parser.add_argument('-m','--mapping-quality',dest='mapq',help="Mapping quality",default=20)
    args = parser.parse_args()
    assert args.bam_file is not None, \
            "Need a bam file as input"
    assert args.reference_name is not None, \
            "Make sure you input a reference, which is the actual referenc rest are contamination spots"
    get_reads(args.bam_file, args.reference_name, mapq=args.mapq)

if __name__=="__main__":
    main()
