#!/bin/bash
# run the best practice gatk analysis for calling variants.

get_options(){
    while getopts "AT:sI:i:pc:mM:r:R:d:mhD" opt; do
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
        A)
            ANCIENT_FASTQ_FILTER="TRUE"
            ;;
        T)
            TRAITS_FILE=$OPTARG  
            ;;
        s)
            END="se"
            ;;
        S)
            SPECIES=$OPTARG
            ;;
        D)
            MAP_DAMAGE="TRUE"
            ;;
        p)
            PMD="TRUE"
            ;;
        m)
            MTDNA="TRUE"
            ;;
        I)
            INDEL_ALIGNMENT=$OPTARG
            ;;
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
        -D use Map Damage.
        -p <PMD THRESHOLD> default: do not use PMD
        -t <TRAITS FILE> default: no traits added to the nexus file
        -h print this output
EOF
}
DIR=$( dirname "$(readlink  $0)")
reference="$DIR/../ref/rCRS.fa"
PICARD="$DIR/../src/picard"
GATK="$DIR/../src/gatk/GenomeAnalysisTK.jar"
RSCRIPTS="$DIR/../rscripts"
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
# Sam post -fix (variable changes throughout the script.
SAM_SEARCH_EXPAND=*.sam
#Read group stuff
END=pe
RGPL=Illumina
get_options "$@"
PATH=$PATH:$DIR/../bin
#gcc -lgfortran
#exit 1
# Default settings if you don't specif anything, 
SETUP_FILE=pipeline_setup.txt
#Specfied the minimum depth 
if [[ ! -f $SETUP_FILE ]]; then
    echo "Cannot find setup file = ${SETUP_FILE}"
    exit 1 
fi
if [[ ! -f .fin_pipeline ]]; then
    rm -Rf temp_results/*
    rm -Rf $results_dir/coverage
fi
mkdir -p temp_results
tmp_dir=temp_results
results_dir=results
#addmap_damage folder
mkdir -p $results_dir
mkdir -p $results_dir/damage
mkdir -p $results_dir/coverage
mkdir -p $results_dir/pmd
echo $SAM_SEARCH_EXPAND
#Source after the environment has been setup
if [[ $PMD != "" ]]; then
    pmd
    index_bams
    exit 1
fi
source $DIR/ancient_dna_funcs.sh
#source $DIR/modern_human_funcs.sh
# This code generates the reference genome that is 2 x the original.
#if [[ $MTDNA != "" ]]; then
#    touch ref.tmp
#    cat $reference > ref.tmp
#    cat $reference | grep -V "^>"  > r.tmp
#    cat r.tmp >> ref.tmp
#    rm r.tmp
#    reference=ref.tmp
#fi
if [[ $ANCIENT_FASTQ_FILTER = "TRUE" ]]; then
    ancient_filter
fi

map_reads 
sort_bam
mark_duplicates
index_bams
#Run some map Damage
#consensus_sequences

# TODO COMPARE HaplotypeCaller and Samtools
#call_variants_samtools
#if [[ $MAP_DAMAGE != "" ]]; then
#    map_damage   
#    index_bams
#fi
add_and_or_replace_groups 
index_bams
if [[ $PMD != "" ]]; then
    pmd
    index_bams
fi
#if [[ $MINIMAL != "TRUE" ]]; then
#    haplotype_caller
#fi
haplocaller_combine
vcf_filter
vcf_to_haplogrep
coverage_plots_R
remove_g_a_c_t
# Turn them all the fasta
vcf_to_fasta
# Post-mortem damage 
##fasta_to_nexus
align_muscle
if [[ $TRAITS_FILE != "" ]]; then
    annotate_traits
fi
touch .fin_pipeline
