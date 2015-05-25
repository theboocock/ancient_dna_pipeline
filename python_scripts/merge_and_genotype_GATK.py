#!/usr/bin/env python
#
# The script is used to mergetheGVCFs into a smaller set
# of file when too many files are present.
#
#

### Boiler plate to get parallel jobs running in python

import argparse
from job_utilities import *

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

SPLIT_SIZE=100
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

def haplotype_caller(gatk, xmx, cores, reference, inputs,output):
   command = HAPLOTYPE_CALLER_TEMPLATE.format(xmx, gatk,  reference, output) 
   command = command + ' --variant ' +' --variant '.join(inputs[0:100])
   queue_jobs([command],'haplotypeCaller',1)

def haplotype_single(gatk, xmx, cores, reference, inputs):
   commands = []
   for sample in inputs:
       output = sample + 'test.vcf'
       command = HAPLOTYPE_CALLER_TEMPLATE.format(xmx, gatk,  reference, output) 
       command = command + ' --variant ' +  sample
       commands.append(command)
   queue_jobs(commands,'haplotypeCaller',1)

def main():
    parser = argparse.ArgumentParser(description='MergeGVCFs and genotype them using the GATK')
    parser.add_argument('-g','--gatk',dest='gatk',help="Location of the GATK", required=True)
    parser.add_argument('-x','--xmx',dest='xmx', help="Memory to use with JAVA", required=True)
    parser.add_argument('-c','--cores', dest='cores', help="Number of cores to use")
    parser.add_argument('-o','--output',dest='output', help='Final output from the haplotype caller')
    parser.add_argument('-r','--reference', dest='reference', help='Reference FASTA file')
    parser.add_argument('gvcfs',nargs="*",help='gVCF variant call files output from the GATK')
    args = parser.parse_args()
    args.cores = int(args.cores)
    # Remove quotes needed to pass the Xmx argument through into the script
    args.xmx= args.xmx.strip('"')
    #outputs = merge_gvcfs(gatk=args.gatk, xmx=args.xmx, cores=args.cores, gvcfs=args.gvcfs, reference=args.reference) 
    haplotype_caller(gatk=args.gatk, xmx=args.xmx, cores=args.cores, inputs=args.gvcfs, output=args.output, reference=args.reference)
    #haplotype_single(gatk=args.gatk, xmx=args.xmx, cores=args.cores, inputs=args.gvcfs, reference=args.reference)


if __name__ == "__main__":
    main()

