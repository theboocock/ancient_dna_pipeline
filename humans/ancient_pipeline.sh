#!/bin/bash
# run the best practice gatk analysis for calling variants.

get_options(){
    while getopts "s:i:c:m:M:r:R:d:mh" opt; do
        case $opt in
        i)
            SETUP_FILE=$OPTARG
            ;;
        m)
            XMX=-Xmx${OPTARG}
            ;;
        M)
            MAPPER=$OPTARG
            ;;
        r) 
            reference=$OPTARG
            ;;
        R)
            RGPL=$OPTARG
            ;;
        d)
            results_dir=$OPTARG
            ;;
        s)
            END="se"
            ;;
        S)
            SPECIES=$OPTARG
        h)
            usage
            exit 1
            ;;
        ?)
            usage
            exit 1
        ;;            
        esac
    done
}

usage(){
cat << EOF
        Usage: Runs DNA processing for the anatomy laboratory.

        -i <Setup File> 
        -c <Cores to use>
        -m <Memory for Java VMs>
        -M <Mapper bwa or bowtie>
        -r <Reference Genome>
        -R <Read platform>
        -d <Output Directory>
        -s <Single ended> default: paired end.
        -S <Species> default: human
        -h print this output
EOF
}
PATH=$PATH:$DIR/../bin
DIR=$( dirname "$(readlink  $0)")
# Default settings if you don't specif anything, 
SETUP_FILE=pipeline_setup.txt
#Specfied the minimum depth 
MIN_DEPTH=10
#Specify the number of cores to use
CORES=6
XMX=-Xmx2g
SPECIES=human
# JAVA7
# mapper of choice either bwa of bowtie at the moment
MAPPER=bwa
JAVA7=/Library/Java/JavaVirtualMachines/jdk1.7.0_60.jdk/Contents/Home/bin/java
# bowtie index requires the basename of the file
# Defaults to the human reference Genome. 
reference="$DIR/../ref/rCRS.fa"
PICARD="$DIR/../src/picard"
GATK="$DIR/../src/gatk/GenomeAnalysisTK.jar"
RSCRIPTS="$DIR/../rscripts"
# Sam post -fix (variable changes throughout the script.
SAM_SEARCH_EXPAND=*.sam
#Read group stuff
END=pe
RGPL=Illumina
mkdir -p temp_results
tmp_dir=temp_results
results_dir=results
#addmap_damage folder
get_options "$@"
mkdir -p $results_dir
mkdir -p $results_dir/damage
mkdir -p $results_dir/coverage
echo $SAM_SEARCH_EXPAND
#Source after the environment has been setup

source $DIR/ancient_dna_funcs.sh
map_reads 
add_and_or_replace_groups 
sort_bam
mark_duplicates 
index_bams
#Run some map Damage
consensus_sequences

# TODO COMPARE HaplotypeCaller and Samtools
#call_variants_samtools

if [[$MINIMAL != "TRUE" ]]; then
    haplotype_caller
fi

haplocaller_combine
vcf_filter
vcf_to_haplogrep
coverage_plots_R

# Turn them all the fasta
if [[ $VCF_TO_FASTA != "" ]]; then
    vcf_to_fasta
fi

if [[ $MAP_DAMAGE != "" ]]; then
    map_damage    
fi
