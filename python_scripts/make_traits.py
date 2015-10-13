#!/usr/bin/env python
#| 
#| This python script takes a traits files, which looks like the following
#|
#| samplename trait1 trait2 etc
#|
#| By finding the traits and breaking them into bit strings, so they can be loaded
#| into the nexus output file
#|

import os
import numpy
import argparse
import sys
import uuid

TRAITS_TEMPLATE="""BEGIN TRAITS;
                   Dimensions NTRAITS={0};
                   Format labels=yes missing=? separator=Comma;
                   TraitLabels {1};
                   Matrix
                   {2}
                   ;
                   END;
                   """

TAXA_TEMPLATE="""#NEXUS
BEGIN TAXA;
DIMENSIONS NTAX={0};
TAXLABELS
{1}
;
END;
"""

def generate_traits_dictionary(trait_file):
    """
        Generate a dictionary of traits
    """
    trait_list=[]
    with open(trait_file) as traits:
        for trait in traits:
            trait_inner_l=[]
            t_line=trait.split('\t')
            for i in range(0,len(t_line)):
                trait_inner_l.append(t_line[i])
            trait_list.append(trait_inner_l)
    return(numpy.array(trait_list))

def comma_sep_binary_string(sample_names,full_sample_list,individual_trait_list):
    """
        Turns into a comma seperated trait file for nexus
    """
    string_dict={}
    trait_dict={}
    individual_trait_list=[s.strip() for s in individual_trait_list]

    individual_trait_list=["UNKNOWN" if s == "" else s for i, s in enumerate(individual_trait_list)]
    print full_sample_list
    for sample, trait in zip(full_sample_list, individual_trait_list):
        trait_dict[sample] = trait 
    uniq_list = set(individual_trait_list)
    temp_string = "0,"*len(uniq_list)
    temp_string =  temp_string[:len(temp_string)-1]
    n_traits = len(uniq_list)
    trait_labels=[]
    for i, item in enumerate(uniq_list):
        trait_labels.append(item)
        mod_string = list(temp_string)
        mod_string[(i)*2] = "1"
        string_dict[item] = ''.join(mod_string)
    sample_dict = {}
    for sample in sample_names:
        sample_dict[sample] = string_dict[trait_dict[sample]]
    return(sample_dict, n_traits,trait_labels)

def traits_to_nexus(input_file ,output_prefix, traits):

    """
        Reads traits and processes the pandas trait dictionary.
    """
    # Pull out the sample names
    samplenames = traits[1:,0]
    temp_file = str(uuid.uuid4())
    with open(temp_file,'w') as out_file:
        a = open(input_file)
        temp_middle=a.read()
        f_ind=temp_middle.find("matrix") + 7
        semi_ind=temp_middle.find(';',f_ind)
        actual_sample_names=[]
        for samps in temp_middle[f_ind:semi_ind].split('\n'):
            try:
                actual_sample_names.append(samps.split()[0])
            except IndexError:
                break
        out_file.write(TAXA_TEMPLATE.format(len(actual_sample_names),"\n".join(actual_sample_names)))
        temp_middle=temp_middle.replace("#NEXUS","")
        out_file.write(temp_middle)
        a.close()
    with open(temp_file,'r') as template_file:
        sequence_and_taxa=template_file.read()
        for i in range(1,traits.shape[1]):
            temp_dict={}
            trait_name = traits[0,i].replace(" ", "")
            t_list = traits[1:,i]
            out_file=open(output_prefix + '.' + trait_name.strip() + '.nex','w')
            out_file.write(sequence_and_taxa)
            (sample_dict, n_traits, trait_labels)=comma_sep_binary_string(actual_sample_names, samplenames,t_list)
            matrix = ""
            for key, value in sample_dict.items():
                matrix += key + " " + value + "\n"
            trait_labels=[o.replace(" ","") for o in trait_labels]
            trait_labels=' '.join(trait_labels)
            out_file.write(TRAITS_TEMPLATE.format(n_traits,trait_labels ,matrix))
    os.remove(temp_file)

def main():
    parser = argparse.ArgumentParser(description="Takes a trait file and creates traits for all the descriptions")
    parser.add_argument('-i','--nex',dest='nexus_input',help='Nexus Input File')
    parser.add_argument('-o','--out-nex',dest='nexus_output',help='Nexus Output File prefix')
    parser.add_argument('-t','--traits',dest='traits_file',help='Traits input File')
    args = parser.parse_args()
    assert args.nexus_input is not None, " A nexus input file is needed for the option --nex or --i"
    assert args.nexus_output is not None, " A nexus output file is needed for the option -o or --output-nex"
    assert args.traits_file is not None, " A nexus traits file is needed for the option -t or --traits"
    try:
        a = open(args.nexus_input)
    except IOError:
        sys.stderr.write('Cannot open nexus_input file: {0}'.format(a))
    try:
        a= open(args.traits_file)
    except IOError:
        sys.stderr.write('Cannot open nexus_output file: {0}'.format(a))

    trait_list=generate_traits_dictionary(args.traits_file)
    traits_to_nexus(args.nexus_input, args.nexus_output, trait_list)

if __name__=="__main__":
    main()
