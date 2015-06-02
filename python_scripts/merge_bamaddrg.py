#!/usr/bin/env python
#
#
# Prepare to run bamdrg



import argparse
import pysam
import subprocess


def create_cmd_line(bams, output):
   cmd = []
   
   cmd.append('bamaddrg')

   for item in bams:
       bamfile = pysam.AlignmentFile(item)
       sample_name = (bamfile.header['RG'][0]['SM'])
       cmd.extend(['-b', item ,'-s',sample_name])
   subprocess.check_call(cmd,stdout=open(output,'w')) 


def main():
    parser = argparse.ArgumentParser(description="Processes bams through bamaddrg")
    parser.add_argument('bams',nargs="*", help="bam files to merge using bamarg")
    parser.add_argument('-o','--output',dest='output_bam',help='Bam file output')
    args = parser.parse_args()
    assert args.output_bam is not None, \
            "Ensure you have a -o or --output argument"
    create_cmd_line(args.bams, args.output_bam)

if __name__=="__main__":
    main()
