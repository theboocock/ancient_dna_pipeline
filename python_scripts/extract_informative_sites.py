#!/usr/bin/env python
#
# Extract informative sites from a VCF
#
# @author James Boocock
# 
#


import pyfasta
import argparse
import sys

def check_fasta(fasta):
    sequences = sorted(fasta.keys())
    # make sure the lengths are the same
    for i, sequence in enumerate(sequences):
        if i == 0:
            length = len(fasta_r[sequence])
            prev_length = len(fasta_r[sequence])
        else:
            length = len(fasta_r[sequence])
            if length != prev_length:
                sys.stderr.write("Fasta sequences are not the same length! Please check fasta") 
                sys.exit(1)


def extract_informative(fasta, indels=None):
    """
        Extract variable sites from a fasta file

        Writes the positions of the variable sites into a file. 
        Currently this only works with imputted fasta files that do not include indels.
        TODO

        Make work with indels maybe extract from a VCF file. Dunno yet 
    """
    sequences = sorted(fasta.keys())
    seq_len = sequences[0]
    
    out_sequences = {}
    for sequence in sequences:
        out_sequences[sequence] = ""
    with open("informative_positions.txt", "w") as f:
        for i in range(len(fasta[seq_len])):
            alleles = []
            skip_allele = False
            for sequence in  sequences:
                tmp_allele = fasta[sequence][i]
                if tmp_allele == "-" or tmp_allele == "N":
                    skip_allele = True
                    break
                alleles += fasta[sequence][i]
                # Ignore indels for now
            if not skip_allele and len(set(alleles)) > 1:
                f.write(str(i + 1) + "\n")
                for j, sequence in enumerate(sequences):
                    out_sequences[sequence] += alleles[j]
    for sequence in sequences: 
        print ">"+sequence
        if indels == None:
            print out_sequences[sequence]
        

def main():
    parser = argparse.ArgumentParser(description="Extract informative sites")
    parser.add_argument("aligned_fasta")
    parser.add_argument("-v","--vcf", dest="vcf", help="Extract genotypes from VCF")
    args = parser.parse_args()
    fasta_r = pyfasta.Fasta(args.aligned_fasta)
    indels=None
    if args.vcf is not None:
        # We should generate a list of indels for each sample
        indels = get_indels_from_vcf(args.vcf)
    extract_informative(fasta_r,indels)

if __name__=="__main__":
    main()

