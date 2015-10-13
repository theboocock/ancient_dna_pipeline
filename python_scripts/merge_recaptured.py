#!/usr/bin/env python
#
#
#
#

import argparse
from subprocess import check_call

def process_fastqs(fastqs):
    fastq_dict = {}
    for fastq in fastqs:
        try:
            name = fastq.split('_')[0]
            name = name.split('-')[0]
            read = fastq.split('_')[3]
            name = name+"_"+read
            try:
                fastq_dict[name].append(fastq)
            except KeyError:
                fastq_dict[name]= [fastq]
        except IndexError:
            pass
    for fq_name, fq_list in fastq_dict.items():
            name = fq_name.split('_')[0]
            read = fq_name.split('_')[1]
            output_name = open(name+ "_merged_libraries_"+read+"_001.fastq", 'w')
            cmd = ['cat']
            cmd.extend(fq_list)
            check_call(cmd, stdout=output_name)
def main():
    parser = argparse.ArgumentParser(description="Merges reads")
    parser.add_argument("fastqs",nargs="*")
    args = parser.parse_args()
    process_fastqs(args.fastqs)

if __name__=="__main__":
    main()
