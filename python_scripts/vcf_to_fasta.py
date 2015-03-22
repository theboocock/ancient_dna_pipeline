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


def vcf_to_fasta(input_vcf, output_fasta, ref_seq, species, use_indels, min_depth,min_probs=0.9):
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
            is_beagle = False
            try:
                pl = sample['PL']
                pheno_l = [int(o) for o in pl]
                dp = sample['DP']
                pl = pheno_l.index(min(pheno_l))
                if(genotype == None or float(dp) <= min_depth):
                    sample_fasta[sample.sample][temp_position] = 'N'
                    # Just to ensure, the bad thing doesn't occur
                    # Overwriting the N call.
                    continue
            except AttributeError:
                is_beagle = True
                gp = sample['GP']
                g_l = [float(o) for o in gp]
                if( max(g_l) < min_probs):
                    sample_fasta[sample.sample][temp_position] = 'N'
                    continue
                pl = g_l.index(max(g_l))
            except TypeError:
                sample_fasta[sample.sample][temp_position] = 'N'
                continue
            position = record.POS
            temp_position = position - 1 
            sample = sample.sample
            if not is_beagle:
                genotype=genotype.split('/')
            else:
                genotype=genotype.split("|")

            # If pl is greater than zero
            ref=record.REF
            alt=record.ALT
            # Gl is substituted
            if(int(pl) >0 ):
                if (is_ga_or_ct(ref, alt)):
                    if(is_beagle):
                        if(g_l[0] > g_l[2]):
                            continue
                    elif(pheno_l[0] < pheno_l[2]):
                        continue
                no_alleles = 1 + len(alt)
                genotype = genotype[0]
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
                        else:
                            gt = real_gt[i]
                            sample_fasta[sample][temp_position] = gt[0]
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
    parser.add_argument('--min-depth',dest="min_depth", 
                        default=5)
    args = parser.parse_args()
    assert  args.fasta_output is not None, \
            "-o or --output is required"
    assert args.vcf_input is not None, \
            "-i or --vcf is required"
    assert args.reference is not None, \
            "-r or --reference is required"
    if(args.use_indels is None):
        args.use_indels = False
    vcf_to_fasta(args.vcf_input, args.fasta_output, args.reference, args.species,args.use_indels, args.min_depth) 
if __name__ == "__main__":
    main()
