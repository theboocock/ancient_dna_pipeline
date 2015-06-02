#!/usr/bin/env python
"""
    Script parser the output of coverage_plots R

"""

import argparse
import subprocess
import shlex
import sys
def _get_samples(coverage_file):
    """
        Private function that returns a list of samples from the coverage data file.
    """
    cov_hash = {}
    with open(coverage_file) as cov:
        for line in cov:
            line = line.split()
            sample = line[0].split('.')[0]
            coverage_percent = float(line[1])
            cov_hash[sample] = coverage_percent
    return cov_hash

def filter_sample(min_coverage, coverage_file):
    """
        Function filters a sample for minimum coverage using a coverage file.
    """
    cov_hash = _get_samples(coverage_file)
    samples_to_keep = []
    for sample, coverage in cov_hash.items():
        if coverage >= min_coverage:
            samples_to_keep.append(sample)
    return samples_to_keep

def bcftools_filter(vcf_input, vcf_output, samples_to_keep):
    """
        Filter vcf using bcftools.
    """
    command = 'bcftools view -s '
    sample_comma = ','.join(samples_to_keep)
    command = command + sample_comma + ' -o ' + vcf_output + ' '
    command += vcf_input + ' '
    try:
        command = shlex.split(command)
        print command
        subprocess.check_call(command)
    except subprocess.CalledProcessError:
        sys.stderr.write('Failed to run bcftools\n')
        sys.exit(10)

def main():
    """
        main function for filter samples file

    """
    #TODO Add lot's more metrics, for now just filter samples from the VCF using coverage.
    #     This should mean we can do a final grep of all the samples.
    parser = argparse.ArgumentParser(description="Read coverage plots file")
    parser.add_argument('-m', '--min-coverage', dest='minimum_coverage',
                        help="Minimum coverage before filtering a sample", default=95)
    parser.add_argument('-c', '--coverage_file', dest='coverage_file',
                        help="Coverage plots", required=True)
    parser.add_argument('-v', '--vcf', dest='vcf_input', help='Vcf input file', required=True)
    parser.add_argument('-o', '--output-vcf', dest='vcf_output',
                        help='Vcf output file', required=True)
    args = parser.parse_args()

    args.minimum_coverage = float(args.minimum_coverage)
    if args.minimum_coverage < 1:
        args.minimum_coverage *= 100

    samples_to_keep = filter_sample(min_coverage=args.minimum_coverage,
                                    coverage_file=args.coverage_file)
    bcftools_filter(vcf_input=args.vcf_input, vcf_output=args.vcf_output,
                    samples_to_keep=samples_to_keep)

if __name__ == "__main__":
    main()
