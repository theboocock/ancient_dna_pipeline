#!/usr/bin/env python
#
# The script is used to mergetheGVCFs into a smaller set
# of file when too many files are present.
#
#

### Boiler plate to get parallel jobs running in python

import argparse
from job_utilities import *

HAPLOTYPE_CALLER="""
    java -jar {0} \
    {1} \
    -T HaplotypeCaller \
    -R {2} \
    -o {3} \
"""

MERGE_GVCFS_TEMPLATE="""
    java -jar {0} \
    {1} \
    -T CombineGVCFs \
    -R {2} \
    -o {3} \
"""

HAPLOTYPE_CALLER_TEMPLATE="""
    java -Djava.io.tmpdir=tmpdir -jar {0} \
    {1} \
    -T GenotypeGVCFs \
    -R {2} \
    -o {3} \
"""

def get_bam_pairs(bams):
    """
        Function returns pairs of bam files because sample ID relies
        on samples being encoded with a '.'

        @bams  A List of Bam files.
    """
    ### TODO Use the read group to merge the reads, far smarter!
    bam_list = {}
    for bam in bams:
        sid = bam.split('.')[0]
        try:
            bam_list[sid].append(sid)
        except KeyError:
            bam_list[sid] = [sid]
    return bam_list

def haplotype_caller(gatk, xmx, cores, reference, bams):
    """
        Function creates gVCFS using the GATK.

        @param gatk the GATK
        @param xmx the memory for java
        @param cores the number of cores on the machine to use
        @param reference the reference genome for alignments
        @param bams the alignment files to call GATK.

        @return gvcfs, a list of the gVCFS to use.
    """
    gvcfs = []
    bam_pairs = get_bam_pairs(bams)
    for bams in bam_pairs:
        output = bams[0].split('.')+'.gvcf'
        command = HAPLOTYPE_CALLER_TEMPLATE.format(xmx, gatk, reference, output)
        command = ' -I ' + ' -I '.join
    return gvcfs

def merge_gvcfs(gatk, xmx, cores, gvcfs, reference):
    """
        Function merges GVCFs using the GATK.

        @param GATK the gatk directory
        @param xmx the memory for java
        @param cores the number of cores on the machine to use
        @param gvcfs genotypeGVCFs
        @param reference reference genome

        @return outputs, a list of the mergedGVCFS
    """
    SPLIT_SIZE=100
    commands = []
    outputs = []
    stdouts = []
    no_groups = (len(gvcfs)/SPLIT_SIZE) + 1  
    for i in range(0,no_groups):
        output = str(i) + '.gvcf'
        outputs.append(output)
        command = MERGE_GVCFS_TEMPLATE.format(xmx, gatk ,reference, output) 
        command = command + '--variant ' + ' --variant '.join(gvcfs[i:(i*SPLIT_SIZE + SPLIT_SIZE)])
        commands.append(command)
    queue_jobs(commands, "mergeGVCFs", cores) 
    return(outputs)

def combine_gvcfs(gatk, xmx, cores, reference, inputs,output):
   command = HAPLOTYPE_CALLER_TEMPLATE.format(xmx, gatk,  reference, output) 
   command = command + ' --variant ' +' --variant '.join(inputs[0:100])
   queue_jobs([command],'combineGVCF',1)

#
#def haplotype_single(gatk, xmx, cores, reference, inputs):
#   commands = []
#   for sample in inputs:
#       output = sample + 'test.vcf'
#       command = HAPLOTYPE_CALLER_TEMPLATE.format(xmx, gatk,  reference, output) 
#       command = command + ' --variant ' +  sample
#       commands.append(command)
#   queue_jobs(commands,'haplotypeCaller',1)

def main():
    parser = argparse.ArgumentParser(description='MergeGVCFs and genotype them using the GATK')
    parser.add_argument('-g','--gatk',dest='gatk',help="Location of the GATK", required=True)
    parser.add_argument('-x','--xmx',dest='xmx', help="Memory to use with JAVA", required=True)
    parser.add_argument('-c','--cores', dest='cores', help="Number of cores to use")
    parser.add_argument('-o','--output',dest='output', help='Final output from the haplotype caller')
    parser.add_argument('-r','--reference', dest='reference', help='Reference FASTA file')
    parser.add_argument('bams',nargs="*",help='gVCF variant call files output from the GATK')
    args = parser.parse_args()
    args.cores = int(args.cores)
    args.xmx= args.xmx.strip('"')
    genovcfs = haplotype_caller(gatk=args.gatk, xmx=args.xmx, cores=args.cores, bams=args.bams, reference=args.reference) 
    outputs = merge_gvcfs(gatk=args.gatk, xmx=args.xmx, cores=args.cores, gvcfs=genovcfs, reference=args.reference) 
    combine_gvcfs(gatk=args.gatk, xmx=args.xmx, cores=args.cores, inputs=args.gvcfs, output=args.output, reference=args.reference)
    #haplotype_single(gatk=args.gatk, xmx=args.xmx, cores=args.cores, inputs=args.gvcfs, reference=args.reference)


if __name__ == "__main__":
    main()

