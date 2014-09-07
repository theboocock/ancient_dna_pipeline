#!/usr/bin/env python
"""
    This program takes a vcf input and converts the VCF
    file to a fasta file containing all the relevant information
    for further analysis. In humans only one such problem exists
    with an indel being called in a different place on the reference
    sequence.


"""
import os,re,vcf,argparse
from pyfasta import Fasta

def vcf_to_fasta(input_vcf, output_fasta, ref_seq, species, use_indels):
    # First part is to get the fasta sequence then atke each position 
    # and then alter the reference as necessary for each sample.
    # Because everyone will have different SNPs.
    f = Fasta(ref_seq)
    # For now this is only going to work with mtDNA sequences,
    # but plan to extend this in the future to full genome
    # gets the full genomes sequences and currently assumes 
    # that the fasta only contains one sequence.
    full_sequence = list(f[f.keys()[0]])
    sample_fasta = {}
    vcf_reader = vcf.Reader(open(input_vcf,'r'),strict_whitespace=True)
    samples = vcf_reader.samples
    for sample in samples:
        sample_fasta[sample] = full_sequence[:]
    for record in vcf_reader:
        for sample in record.samples:
            genotype = sample['GT']
            pl = sample['PL']
            if(genotype == None):
                continue
            sample = sample.sample
            position = record.POS
            genotype=genotype.split('/')
            pl = [int(o) for o in pl]
            pl = pl.index(min(pl))
            # If pl is greater than zero
            if(int(pl) > 0):
                alt=record.ALT
                no_alleles = 1 + len(alt)
                ref=record.REF
                genotype = genotype[0]
                temp_position = position - 1 
                real_gt =str(alt[int(genotype)-1])
                if(species == 'human'):
                    if(position == 8270 and ref=="CACCCCCTCT"):
                        sample_fasta[sample][8280:8289] = '-'*9 
                        continue
                for i in range(0,max(len(real_gt),len(ref))):
                    if ( i == (len(real_gt) - 1) and i == (len(ref)- 1)):              
                        gt = real_gt[i]
                        sample_fasta[sample][temp_position] = gt
                    elif(len(real_gt) > len(ref) and i != 0):
                        if(use_indels):
                            gt = list(real_gt[i])
                            sample_fasta[sample]= sample_fasta[sample][:temp_position] + gt + sample_fasta[sample][temp_position:] 
                            temp_position = temp_position + 1 
                    elif(len(real_gt) < len(ref) and i != 0):
                        sample_fasta[sample][temp_position + i] = '-'
    with open(output_fasta,'w') as out:
        for sample in samples:
            out.write('>'+sample + '\n')
            out.write("".join(sample_fasta[sample])+'\n')
   
        
def main():
    parser= argparse.ArgumentParser()
    parser.add_argument("-i",'--vcf',dest="vcf_input",
                        help="VCF input file to convert to fasta")
    parser.add_argument('-o','--output',dest="fasta_output",
                        help="Fasta output file")
    parser.add_argument('-r','--reference',dest='reference',
                        help="Reference FASTA sequence file")
    parser.add_argument('-s','--species',dest='species',default='human',
                        help="Species that you are performing analysis on, "
                             "Currently accepted values are human and dog")
    parser.add_argument('--use-indels',dest='use_indels',action="store_true",
                        help="Do not use indels in the analysis", default=False)
    args = parser.parse_args()
    assert  args.fasta_output is not None, \
            "-o or --output is required"
    assert args.vcf_input is not None, \
            "-i or --vcf is required"
    assert args.reference is not None, \
            "-r or --reference is required"
    if(args.use_indels is None):
        args.use_indels = False
    vcf_to_fasta(args.vcf_input, args.fasta_output, args.reference, args.species,args.use_indels) 
if __name__ == "__main__":
    main()
