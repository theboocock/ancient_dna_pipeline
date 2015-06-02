#!/bin/bash
#
# Indexes a FASTA reference file
# $1 REFERENCE FASTA_FILE

DIR=$( dirname $0)
echo $DIR
PICARD="$DIR/../src/picard/"
DICT_OUTPUT=${1%.fa}.dict
bwa index $1
samtools faidx $1
java -jar -Xmx2g $PICARD/CreateSequenceDictionary.jar  R=$1 O=$DICT_OUTPUT
