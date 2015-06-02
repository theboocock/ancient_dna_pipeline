#!/bin/bash
#
# Generate the pipeline file for both the processed files and the raw fastq files.
# Dynamic Trim using either the dynamic trim
#
#
# $1 NZGL folder fastq root folder 
#$2 processed or raw
#
DATA_FOLDER=$1
PROCESSED_OR_TRIMMED=$2

usage(){
cat << EOF 
    Get the raw files from a NZGL run  
    ./create_pipeline_setup.sh <SEQUENCE_FOLDER> <PROCESSED_OR_TRIMMED>  
EOF
}
if [[ -z DATA_FOLDER ]]; then
    usage
    exit 1
fi

if [[ -z $PROCESSED_OR_TRIMMED ]]; then
    usage
    exit 1
fi
# This is the folder of raw FQ sequences that we get from
# NZGL, if we perform a cut on the files in here we can
# get the sample names for further analysis.


if [[ $OUTPUT = "" ]];  then
    OUTPUT=pipeline_setup.txt
fi

#ls $DATA_FOLDER
if [[ $PROCESSED_OR_TRIMMED = 'trimmed' ]]; then
    ls $1/*R1*trimmed | cut -d '_' -f 2 | xargs basename > samples_list.txt
    for i in $1/*_R1_001.fastq.trimmed ;  do echo ${i} `echo ${i} | sed  's/R1_001/R2_001/g'` >> line_setup.txt; done 
elif [[ $PROCESSED_OR_TRIMMED = 'processed' ]]; then
    ls $1/*R1*fastq | cut -d '_' -f 2 |  xargs basename > samples_list.txt
    for i in $1/*_R1_001.fastq ; do echo ${i} `echo ${i} | sed 's/R1_001/R2_001/g'` >> line_setup.txt; done
elif
    [[ $PROCESSED_OR_TRIMMED = 'raw' ]]; then
     ls $1/*R1*fastq | cut -d '_' -f 1 | sed "s/^\.\///g"  > samples_list.txt
     for i in $1/*_R1_001.fastq ; do echo ${i} `echo ${i} | sed 's/R1_001/R2_001/g'` >> line_setup.txt; done
fi
paste -d ' ' samples_list.txt line_setup.txt > $OUTPUT
rm line_setup.txt samples_list.txt


