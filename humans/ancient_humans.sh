#map ancient_human_data.
#
# author
# James Boocock
#
# read standard and file for pe reads.
#
# for mtDNA samples.
# file looks like this.
# <read_name> <read_1> <read_2>
# $1 input_file
# $2 pe or se 
# $3 results_dir

usage(){
cat << EOF
        Usage: ancient_humans <input_file_name> <pe or se> <result_dir>
EOF

}

parallel_middle_bit(){
    parallel -j 4 "samtools sort {} ${result_dir}/{/.}.sorted" ::: ${result_dir}/*.bam
    parallel -j 4 "java -Xmx1gb -jar ${PICARD}/MarkDuplicates.jar INPUT={} OUTPUT=${result_dir}/{/.}.rmdup.bam METRICS_FILE=${result_dir}/{/.}.mark_dups.log AS=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT" ::: $result_dir/*.sorted.bam
}

if [ "$1" == "" ]; then
    echo "You must specify an input file name "
    usage
    exit 1
fi

PATH=$PATH:$DIR/../bin
DIR=$( dirname "$(readlink  $0)")

# bowtie index requires the basename of the file
rCRS="$DIR/../ref/rCRS"
PICARD="$DIR/../src/picard/"
GATK="$DIR/../src/gatk/GenomeAnalysisTK.jar"
#Read group stuff
RGPL=Illumina
if [ "$3" == "" ]; then
    mkdir -p results
    result_dir=results
else
    mkdir -p $3
    result_dir=$3
fi
while read line
do
 #   echo $line
    file_name=$(echo $line | cut -f 1 -d ' ')
    output=$file_name 
  #  echo $PWD
    if [ $2 == "pe" ]; then
        pe_one=`echo ${line} | cut -d ' ' -f 2`
        pe_two=`echo ${line} | cut -d ' ' -f 3`
        bowtie2 -x $rCRS -1 $pe_one -2 $pe_two > $result_dir/$output.sam 2> $result_dir/$output.bowtie2.err 
    else
        se=$(echo $line | cut -f -d ' ' -f 2)
        bowtie2 -x $rCRS -u $se > $result_dir/$output.sam 2> $result_dir/$output.bowtie.err
    fi  
    java -Xmx1gb -jar ${PICARD}/AddOrReplaceReadGroups.jar INPUT=${result_dir}/$output.sam OUTPUT=$result_dir/$output.bam RGPL=$RGPL RGPU=$file_name RGSM=$file_name 
    # samtools stuff
done < $1
parallel_middle_bit


