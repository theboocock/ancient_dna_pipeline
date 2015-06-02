#!/usr/bin/env python
# Recal VCF, using likelihood column.
# Only if it's a C-T or G-A transition.
#
# This tool is used to zero the info in a VCF file, so that we can then convert the file.
#
# @author: James Boocock
# @date: 16/03/2015
#

import argparse
import sys
import vcf
import collections
import copy

def is_ga_or_ct(ref,alt):
    if(len(ref) == 1 and len(alt) == 1):
        alt = alt[0]
        if(ref == "C" and alt == "T"):
            return True
        elif(ref == "T" and alt == "C"):
            return True
        elif(ref == "G" and alt == "A"):
            return True
        elif(ref == "A" and alt == "G"):
            return True
    else:
        return False

def recal_vcf(input_vcf,min_dp=2):
    vcf_reader = vcf.Reader(open(input_vcf,'r'),strict_whitespace=True)
    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
    for record in vcf_reader:
        temp_record = record
        for i, sample in enumerate(record.samples):
            idx = 0 
            ref = record.REF
            alt = record.ALT
            #record.samples[i].data = collections.namedtuple("CallData",f_keys)
           # print(sample)
            f_keys = (record.samples[i].data._fields)
            f_vals = [ record.samples[i][vx] for vx in (f_keys)]
            handy_dict = dict(zip(f_keys,f_vals))
            depth = handy_dict['DP']
            if (depth is not None):
                if(int(depth) > min_dp):
                    new_values = [ handy_dict['GT'] if  x == 'GT' else None for x in f_keys]
                else:
                    new_values = [ None for x in f_keys ] 
            else:
                new_values = [ handy_dict['GT'] if  x == 'GT' else None for x in f_keys]
            record.samples[i].data = record.samples[i].data._make(new_values)
        vcf_writer.write_record(record)

def main():
    parser = argparse.ArgumentParser(description="Zero out info in VCF file")
    parser.add_argument("input_vcf",help="Input VCF that we are going to recalfrom")
    parser.add_argument("-m","--mapping_output",dest="",help="")
    args = parser.parse_args()
    recal_vcf(args.input_vcf)


if __name__=="__main__":
    main()
