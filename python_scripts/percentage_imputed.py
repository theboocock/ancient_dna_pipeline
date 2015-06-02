#!/usr/bin/env python
#
# Takes a VCF file that is to be imputed into and tell you the number of sites and the total number of sites that are required to be imputed in each sampe
# 
#

import vcf
import argparse

def impute_table(vcf_input):
    vcf_reader = vcf.Reader(open(vcf_input,'r'),strict_whitespace=True)
    sample_count={}
    samples = vcf_reader.samples
    snp_number = 0
    for sample in samples:
        sample_count[sample] = 0
    for record in vcf_reader:
        snp_number += 1
        for sample in record.samples:
            genotype = sample['GT']
            if(genotype is not None):
                sample_count[sample.sample] += 1
    for sample, counts in sample_count.items():
        print sample, 1-float(counts)/snp_number
        


def main():
    parser = argparse.ArgumentParser(description="Number of impute sites per sample")
    parser.add_argument("vcf_input")
    args = parser.parse_args()  
    impute_table(args.vcf_input)

if __name__=="__main__":
    main()

