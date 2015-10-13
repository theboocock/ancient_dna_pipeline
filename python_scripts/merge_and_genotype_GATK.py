#!/usr/bin/env python

"""
 The script is used to mergetheGVCFs into a smaller set
 of file when too many files are present.

 Boiler plate to get parallel jobs running in python
"""

import argparse
import pysam
from job_utilities import *

HAPLOTYPE_CALLER = """
    java -jar {0} \
    {1} \
    -T HaplotypeCaller \
    --emitRefConfidence GVCF
    -R {2} \
    -o {3} \
    --sample_ploidy {4}
"""

MERGE_GVCFS_TEMPLATE = """
    java -jar {0} \
    {1} \
    -T CombineGVCFs \
    -R {2} \
    -o {3} \
"""

GENOTYPEGVCFS_TEMPLATE = """
    java -Djava.io.tmpdir=tmpdir -jar {0} \
    {1} \
    -T GenotypeGVCFs \
    -R {2} \
    -o {3} \
    --standard_min_confidence_threshold_for_calling 10 \
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
        sam = pysam.AlignmentFile(bam,'rb')
        sample_id  = (sam.header['RG'][0]['SM'])
        try:
            bam_list[sample_id].append(bam)
        except KeyError:
            bam_list[sample_id] = [bam]
    return bam_list

def haplotype_caller(gatk, xmx, reference, bams, cores, out_directory, ploidy, bed_file=None):
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
    commands = []
    try:
        os.mkdir(out_directory)
    except OSError:
        pass
    for sample, bams in bam_pairs.items():
        output = os.path.join(out_directory, os.path.basename(sample + '.g.vcf'))
        command = HAPLOTYPE_CALLER.format(xmx, gatk, reference, output, ploidy)
        command = command + ' -I ' + ' -I '.join(bams) 
        command = command + ' -bamout ' + output + ".bam"
        if bed_file is not None:
            command  = command + " -L " + bed_file
        commands.append(command)
        print command
        gvcfs.append(output)
    queue_jobs(commands, "haplotypeCaller", cores)
    return gvcfs

SPLIT_SIZE = 100
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
    no_groups = (len(gvcfs)/SPLIT_SIZE) + 1
    for i in range(0, no_groups):
        output = str(i) + '.g.vcf'
        outputs.append(output)
        command = MERGE_GVCFS_TEMPLATE.format(xmx, gatk, reference, output)
        command = command + '--variant ' + ' --variant '.join(gvcfs[i:(i*SPLIT_SIZE + SPLIT_SIZE)])
        commands.append(command)
    queue_jobs(commands, "mergeGVCFs", cores)
    return outputs


#
#def haplotype_single(gatk, xmx, cores, reference, inputs):
#   commands = []
#   for sample in inputs:
#       output = sample + 'test.vcf'
##       command = HAPLOTYPE_CALLER_TEMPLATE.format(xmx, gatk,  reference, output) 
#       command = command + ' --variant ' +  sample
#       commands.append(command)
#   queue_jobs(commands,'haplotypeCaller',1)

def genotype_gvcfs(gatk, xmx, cores,
        inputs, output,
        reference, bed_file=None):
    """
        Genotype GVCFs using the GATK
    """
    commands = []
    command = GENOTYPEGVCFS_TEMPLATE.format(xmx, gatk, reference, output)
    command = command + ' --variant ' + ' --variant '.join(inputs)
    if bed_file is not None:
        command  = command + " -L " + bed_file
    commands.append(command)
    output = os.path.join(os.path.dirname(output), 'all_sites.vcf')
    command = GENOTYPEGVCFS_TEMPLATE.format(xmx, gatk, reference, output)
    command = command + ' --variant ' + ' --variant '.join(inputs)
    command = command + ' --includeNonVariantSites'
    if bed_file is not None:
        command  = command + " -L " + bed_file
    commands.append(command)
    queue_jobs(commands, "genotypeGVCFs", cores)

def main():
    """
        Main module runs the GATK for the analysis.
    """
    parser = argparse.ArgumentParser(description='MergeGVCFs and genotype them using the GATK')
    parser.add_argument('-g', '--gatk', dest='gatk', help="Location of the GATK", required=True)
    parser.add_argument('-x', '--xmx', dest='xmx', help="Memory to use with JAVA", required=True)
    parser.add_argument('-c', '--cores', dest='cores', help="Number of cores to use")
    parser.add_argument('-o', '--output', dest='output', 
                        help='Final output from the haplotype caller')
    parser.add_argument('-r', '--reference', dest='reference', 
                        help='Reference FASTA file')
    parser.add_argument('-b','--bed', dest='bed_file',
                        help="Bed file for limiting the GATK")
    parser.add_argument('-p', '--ploidy', dest='ploidy', 
                        help="Sample ploidy", default=2)
    parser.add_argument('-d', '--out_directory', dest='directory', help='Output director')
    parser.add_argument('bams', nargs="*", help='gVCF variant call files output from the GATK')
    args = parser.parse_args()
    args.cores = int(args.cores)
    args.xmx = args.xmx.strip('"')
    genovcfs = haplotype_caller(gatk=args.gatk, xmx=args.xmx, cores=args.cores,
                                bams=args.bams, reference=args.reference,
                                out_directory=args.directory, ploidy=args.ploidy, bed_file=args.bed_file)
    outputs = merge_gvcfs(gatk=args.gatk, xmx=args.xmx, cores=args.cores,
                          gvcfs=genovcfs, reference=args.reference)
    genotype_gvcfs(gatk=args.gatk, xmx=args.xmx, cores=args.cores,
                   inputs=outputs, output=args.output, reference=args.reference,bed_file=args.bed_file)
    #haplotype_single(gatk=args.gatk, xmx=args.xmx, cores=args.cores,
                  #   inputs=args.gvcfs, reference=args.reference)



if __name__ == "__main__":
    main()

