#!/usr/bin/env python
#
# The following program removes misbehaving SNPs.
#
#
#

import argparse
import os
import sys
import pysam

def get_partners(bam, bams):
    """
        Returns the partner bam files having the same ID
    """
    partners = []
    for bam_iter in bams:
        if bam in bam_iter:
            partners.append(bam_iter)
    return partners

def remove_files(bams, remove=False):
    """
        Takes a list of bams to remove from the analysis and removes them
    """
    if not remove:
        try:
            os.mkdir('removed_bams')
        except OSError:
            pass
        for bam in bams:
            sys.stderr.write(bam + " file is being moved\n")
            try:
                os.rename(bam,'removed_bams/')
            except OSError:
                sys.stderr.write("Could not move bam file %s \n" % bam)
        

def get_coverage(bams, reference, chromosome, merged, length, depth, proportion):
    """
        Uses pysam tools to create mpileup, and determines the total coverage
        for either pairs of merged files, or individual merged files.
    """
    removed_bams = []
    if merged:
       # Must be using pairs of matched files
        for b1, b2 in bams:
            total_bases = 0
            reads_one = pysam.mpileup('-D', b1)
            reads_two = pysam.mpileup('-D', b2)
            line_one = reads_one[0]
            line_two = reads_two[0]
            line_one_i = 1
            line_two_i = 1
            while True: 
                if chromosome  not in line_one:
                    try:
                        line_one = reads_one[line_one_i] 
                    except IndexError:
                        break
                    line_one_i += 1
                    continue
                elif chromosome not in line_two:
                    try:
                        line_two = reads_two[line_two_i]
                    except IndexError:
                        break
                    line_two_i += 1
                    continue
                lone_split = line_one.split('\t')
                p1 = int(lone_split[1])
                ltwo_split = line_two.split('\t')
                p2 = int(ltwo_split[1])
                if p1 > p2:
                    cov = int(ltwo_split[3])
                    if cov >= depth:
                        total_bases += 1
                    try:
                        line_two = reads_two[line_two_i]
                    except IndexError:
                        break
                    line_two_i += 1
                elif p2 < p1:
                    cov = int(lone_split[3])
                    if cov >= depth:
                        total_bases += 1
                    try:
                        line_one = reads_one[line_one_i]
                    except IndexError:
                        break
                    line_one_i += 1
                elif p1 == p2:
                    cov = int(lone_split[3]) + int(ltwo_split[3])
                    if cov >= depth:
                        total_bases += 1
                    try:
                        line_one = reads_one[line_one_i]
                    except IndexError:
                        break
                    line_one_i += 1
                    try:
                        line_two = reads_two[line_two_i]
                    except IndexError:
                        break
                    line_two_i += 1
            while True:
                if chromosome  not in line_one:
                    cov = int(lone_split[3])
                    if cov >= depth:
                        total_bases += 1
                    try:
                        line_one = reads_one[line_one_i]
                    except IndexError:
                        break
                    line_one_i += 1
                    continue
                cov = int(lone_split[3])
                try:
                    line_one = reads_one[line_one_i]
                except IndexError:
                    break
                if cov >= depth:
                    total_bases += 1
                line_one_i += 1

            while line_two:
                if chromosome  not in line_two:
                    try:
                        line_two = reads_two[line_two_i]
                    except IndexError:
                        break
                    line_two_i += 1
                    continue
                try:
                    line_two = reads_two[line_two_i]
                except IndexError:
                    break
                cov = int(ltwo_split[3])
                line_two = reads_two.readline()
                if cov >= depth:
                    total_bases += 1
                line_two_i += 1
            if float(total_bases)/length < proportion:
                print b1, b2
                removed_bams.extend([b1, b2])
            # Might need to implement simple one-pass strategy here
    else:
        for b in bams:
            total_bases = 0
            reads = pysam.mpileup('-D', b)
            for r in reads:
                if chromosome in r:
                    r_cov = int(r.strip().split('\t')[3])
                    if r_cov >= depth:
                        total_bases += 1
            if float(total_bases)/length < proportion:
                print b
                removed_bams.append(b)
    remove_files(removed_bams)

def verify_merged(bams):
    """
        Function will remove bam files that were not processed until this point.
    """
    filtered_bams = []
    for bam in bams:
        if "no_collapse" in bam:
            print(bam)
            bam_split = os.path.basename(bam).split('.')[0]
            partners = get_partners(bam_split, bams)
            if len(partners) < 2:
                print bam_split, "is going to be filtered because partner is missing," \
                        "try re-running pipeline with the -M option set if you beleieve this is in error"
            else:
                filtered_bams.append(partners)
    filtered_bam_check = [item for l in filtered_bams for item in l]
    removed_bams = [bam for bam in bams if bam not in filtered_bam_check]
    remove_files(removed_bams)
    return filtered_bams

def main():
    """
        creates argument parser for pipeline filter.
    """

    parser = argparse.ArgumentParser(description="Filters and can remove samples that are "
                                     "misbehaving from the pipeline after alignment")
    parser.add_argument("bams", nargs='*')
    parser.add_argument("-m", "--merged", dest="merged", action="store_true", default=False,
                        help="Defaults to one read per chromosome utilised")
    parser.add_argument("-r", "--reference", dest="reference", required=True)
    parser.add_argument("-c", "--chromosome", dest="chromosome")
    parser.add_argument("--min-coverage", dest="coverage_minimum", default=0.95,
                        help="The minimum coverage percentage for the chromosome of interest")
    parser.add_argument("-s", "--species", dest="species", default="human",
                        help="Species for the analysis")
    parser.add_argument('-l', '--length', dest="length",
                        default="length", help="Length of chromosomal read")
    parser.add_argument('-d', '--depth', dest="depth",
                        default=1, help="Minimum depth to count it being covered by a base")
    parser.add_argument('-p', '--proportion', dest="proportion",
                        default=0.95, help="Samples with coverage less than this are removed")
    ## TODO Add filtering on other metrics, such as Mapping quality
    args = parser.parse_args()
    if args.merged:
        args.bams = verify_merged(args.bams)
    if args.chromosome is None:
        args.chromosome = "gi|251831106|ref|NC_012920.1|"
        args.length = 16569
        if args.species != 'human':
            sys.stderr.write("Cannot default chromosome on anything but human mtDNA\n")
            sys.exit(1)
    args.length = int(args.length)
    args.depth = int(args.depth)
    args.proportion = float(args.proportion)
    get_coverage(args.bams, args.reference, args.chromosome, args.merged, args.length, args.depth, args.proportion) 

if __name__ == "__main__":
    main()
