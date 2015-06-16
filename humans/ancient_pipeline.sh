#!/bin/bash
# run the best practice gatk analysis for calling variants.

get_options(){
    while getopts "tC:AT:sI:i:pc:mMr:R:d:mhDb:P:S:" opt; do
        case $opt in
        P)
            echo $OPTARG
            PLOIDY=$OPTARG
            ;;
        t)
            IMPUTATION="TRUE"
            ;;    
        b)
            BAIL_POS=$OPTARG
            ;;
        i)
            SETUP_FILE=$OPTARG
            ;;
        c)
            CORES=$OPTARG
            ;;
        m)
            XMX=-Xmx${OPTARG}
            ;;
        M)
            MERGED_READS_ONLY="TRUE"
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
            START_POS=$OPTARG
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
        C)
            CONTAMINATION_MAPPING=$OPTARG
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
        -C Using contamination mapping, this is the sequence to extract
        -P <Ploidy> : ploidy value to use for the haplotype caller.
EOF
}

get_params(){
    # AdapterRemoval
    AD_MIN_LENGTH=25
    AD_OVERLAP_LENGTH=11
    # HaplotypeFilter
    HC_QD='2'
    HC_FS='60'
    HC_MQ='30'
    HC_MQRankSum='-12.5'
    HC_ReadPosRankSum='-8.0'
    # BWA
    BWA_ANCIENT_ARG="-n 0.03 -o 2 -l 1024"
}


# This lets us ignore the 

DIR=$( dirname "$(readlink  $0)")
reference="$DIR/../ref/rCRS.fa"
PICARD="$DIR/../src/picard"
GATK="$DIR/../src/gatk/GenomeAnalysisTK.jar"
BEAGLE="$DIR/../src/beagle/beagle.jar"
RSCRIPTS="$DIR/../rscripts"
MAX_FILES=100
MIN_DEPTH=2
#Specify the number of cores to use
CORES=6
XMX=-Xmx16g
# JAVA7
# mapper of choice either bwa of bowtie at the moment
MAPPER=bwa
JAVA7=java
# bowtie index requires the basename of the file
# Defaults to the human reference Genome. 
# Sam post -fix (variable changes throughout the script.
SAM_SEARCH_EXPAND=*.sam
#Read group stuff
END=pe
RGPL=Illumina
#TEST_FREEBAYES="TRUE"
#START_POS="SKIP_READS"
START_POS='MAP_READS'
PLOIDY=1
get_options "$@"
PATH=$PATH:$DIR/../bin
#gcc -lgfortra
echo $CONTAMINATION_MAPPING
echo $SPECIES
if [[ $SPECIES = "" ]]; then 
    echo "You need to specify a species, valid values are human and dog"
    exit 1
fi
if [[ $TEST_FREEBAYES = "TRUE" ]]; then
    echo "We have testing freebayse working"
fi

if [[ $CONTAMINATION_MAPPING != "" ]]; then
    echo "We have contamination mapping working"
    echo $CONTAMINATION_MAPPING  
fi

if [[ $IMPUTATION = "TRUE" ]]; then
    echo "Imputation is working"
fi 
if [[ $MERGED_READS_ONLY = "" ]]; then
    echo "We are not just using the merged reads"
else
    echo "Using only the merged reads"
fi
echo $MAP_DAMAGE
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
mkdir -p $results_dir/bams
mkdir -p $results_dir/contamination
mkdir -p $results_dir/merged
mkdir -p $results_dir/map_damage_plots
echo $SAM_SEARCH_EXPAND
#Source after the environment has been setup
if [[ $PMD != "" ]]; then
    pmd
    index_bams
    exit 1
fi
source $DIR/ancient_dna_funcs.sh
source $DIR/fastq_filters.sh
#source $DIR/modern_human_funcs.sh
# This code generates the reference genome that is 2 x the original.
#if [[ $MTDNA != "" ]]; then
#    touch ref.tmp
#    cat $reference > ref.tmp
#    cat $reference | grep -V "^>"  > r.tmp
#    cat r.tmp >> ref.tmp
#    rm r.tmp
#    reference=ref.tmp
touch .fin_pipeline
if [[ $ANCIENT_FASTQ_FILTER = "TRUE" ]]; then
    ancient_filter
fi

echo "MApdamage status"
echo $MAP_DAMAGE

if [[ $START_POS = 'MAP_READS' ]]; then
    map_reads
    echo "DONE MAP READS" >> .fin_pipeline
    sort_bam
    echo "DONE SORT BAM" >> .fin_pipeline
    if [[ $MAP_DAMAGE != "TRUE" ]]; then
        mark_duplicates
        echo "DONE MARK DUPLICATES" >> .fin_pipeline
    fi
    index_bams
    add_and_or_replace_groups 
    echo "DONE REPLACE_GROUPS" >> .fin_pipeline
    index_bams
    echo "DONE INDEX BAMS" >> .fin_pipeline
    if [[ $CONTAMINATION_MAPPING != "" ]]; then
        save_contaminants
        remove_contaminants
    fi
    store_bams
    echo "DONE STORE BAMS" >> .fin_pipeline
    index_bams
fi



SAM_SEARCH_EXPAND="${results_dir}/bams/*.bam"
merge_bams
#remove_bad_samples
#merge_the_same_samples

# TODO up until this point we are running single samples individually. 
# Files may end up missing in the bams folder so need to be taken care
# of before running into the next step
#Run some map Damage
# TODO COMPARE HaplotypeCaller and Samtools
#call_variants_samtools
if [[ $MAP_DAMAGE != "TRUE" ]]; then
    map_damage  
    echo "DONE MAP DAMAGE" >> .fin_pipeline 
    index_bams
    echo "DONE INDEX BAMS" >> .fin_pipeline
fi
if [[ $PMD != "" ]]; then
    pmd
    echo "DONE PMD" >> .fin_pipeline
    index_bams
    echo "DONE INDEX BAMS" >> .fin_pipeline
fi
if [[ $MINIMAL = "TRUE" ]]; then
    haplotype_caller
    echo "DONE HAPLOTYPECALLER" >>.fin_pipeline
fi
haplocaller_combine
echo "DONE HAPLOCALLER COMBINE" >> .fin_pipeline

if [[ $MAP_DAMAGE != "" ]]; then
    contamination_percentage
    echo "DONE COVERAGE_PLOTS" >> .fin_pipeline
fi
coverage_plots_R

if [[ $IMPUTATION = "TRUE" ]]; then
    # Imputation consists of two distinct steps,
    # Recalling the VCF, then using that with beagle imputation
    #
    beagle_imputation
fi


vcf_filter
echo "DONE VCF FILTER" >> .fin_pipeline
vcf_to_snp_list

vcf_to_haplogrep
echo "DONE VCF HAPLOGREP" >> .fin_pipeline

vcf_to_fasta
# VCF_to_fasta_before muscle
align_muscle
# Post-mortem damage 
fasta_to_nexus
if [[ $TRAITS_FILE != "" ]]; then
    annotate_traits
fi

map_damage_filtered_plots
 Clear this fucknig tmp_dir

CLEAR_DIR="TRUE"
if [[ $CLEAR_DIR = "TRUE" ]]; then 
    rm -Rf $tmp_dir
fi
