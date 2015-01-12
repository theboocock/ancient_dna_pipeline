#!/usr/bin/env python
#
# Covert VCF to R for SNP analysis
#
# POS REF ALT SAMPLE 1 ... SAMPLEN
# 1   A   T   0/0(ref)   ....      1/1(alt)
#
#


import argparse

def get_genotype(gt):
    """ Get the genotypes
        0/0 = 0
        0/1 = 1
        1/1 = 2
        ? What happens if there a non-biallelic SNPs
    """
    if ('0/0'):
        return "0"
    elif ('0/1'):
        return "1"
    else:
        return "2"

def vcf_to_r(vcf_input):
    with open(vcf_input) as vcf_i:
        for line in vcf_i:
            if "##" in line:
                continue
            else:
                #Get header
                if "#" in line:
                    line = line.strip().split('\t')
                    pos = line[1]
                    id = line[2]
                    ref = line[3]
                    alt = line[4]
                    samples = line[9:]
                    print('\t'.join([pos,id,ref,alt]+samples))
                else: 
                    line = line.strip().split('\t')
                    pos = line[1]
                    id = line[2]
                    ref = line[3]
                    alt = line[4]
                    gts = line[9:]
                    gts = [get_genotype(o.split(':')[0]) for o in gts]
                    print('\t'.join([pos,id,ref,alt]+gts))


def main():
    parser = argparse.ArgumentParser(description="VCF to R")
    parser.add_argument('input_file')
    args = parser.parse_args()
    vcf_to_r(args.input_file)

if __name__=="__main__":
    main()
