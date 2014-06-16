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
#
usage(){
cat << EOF
        Usage: ancient_humans <input_file_name> <pe or se> <tmp_dir>
EOF

}

sort_bam(){    
parallel -j ${CORES} "samtools sort {} ${tmp_dir}/{/.}.sorted" ::: ${tmp_dir}/*.bam
    #parallel -j 6 "java -Xmx1g -jar ${PICARD}/ValidateSamFile.jar INPUT={} OUTPUT=${tmp_dir}/{/.}.final.txt" ::: ${tmp_dir}/*.rmdup.bam
}
mark_duplicates(){
    parallel -j ${CORES} "java -Xmx1g -jar ${PICARD}/MarkDuplicates.jar INPUT={} OUTPUT=${tmp_dir}/{/.}.rmdup.bam METRICS_FILE=${tmp_dir}/{/.}.mark_dups.log AS=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT" ::: $tmp_dir/*.sorted.bam
}
consensus_sequences(){
    parallel -j ${CORES} "samtools mpileup -uf ${rCRS}.fa {} | bcftools view -cg - | vcfutils.pl vcf2fq | sed '/^@/!d;s//>/;N' > ${results_dir}/{/.}.fa" ::: ${tmp_dir}/*.rmdup.bam 
}
call_variants(){
    vcf_output=$run_setup_file.vcf
    samtools mpileup -uf ${rCRS}.fa ${tmp_dir}/*.rmdup.bam  | bcftools view -cg - | vcfutils.pl varFilter -d ${MIN_DEPTH}> $results_dir/$vcf_output 
}
vcf_to_haplogrep(){
    vcf_to_haplogrep.py -i $results_dir/$vcf_output -o $results_dir/haplogrep.txt
}

#make_coverage_plots.py has hardcoded settings for coverage

coverage_plots(){
    parallel -j ${CORES} "samtools mpileup -D {} > $tmp_dir/{/.}.cov" ::: ${tmp_dir}/*.rmdup.bam
    make_coverage_plots.pl -- $tmp_dir/*.cov > $results_dir/coverage_plots.ps
}
if [ "$1" == "" ]; then
    echo "You must specify an input file name "
    usage
    exit 1
fi

PATH=$PATH:$DIR/../bin
DIR=$( dirname "$(readlink  $0)")

#input setup file 
run_setup_file=$1

#Specfied the minimum depth 
MIN_DEPTH=10

#Specify the number of cores to use
CORES=6


# bowtie index requires the basename of the file
rCRS="$DIR/../ref/rCRS"
PICARD="$DIR/../src/picard/"
GATK="$DIR/../src/gatk/GenomeAnalysisTK.jar"
#Read group stuff
RGPL=Illumina
mkdir -p temp_results
tmp_dir=temp_results
if [ "$3" == "" ]; then
    mkdir -p results
    results_dir=results
else
    mkdir -p $3
    results_dir=$3
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
        bowtie2 -x $rCRS -p 6 -1 $pe_one -2 $pe_two > $tmp_dir/$output.sam 2> $tmp_dir/$output.bowtie2.err  
    else
        se=$(echo $line | cut -f -d ' ' -f 2)
        bowtie2 -x $rCRS -u $se > $tmp_dir/$output.sam 2> $tmp_dir/$output.bowtie.err
    fi  
    java -Xmx1g -jar ${PICARD}/AddOrReplaceReadGroups.jar INPUT=${tmp_dir}/$output.sam OUTPUT=$tmp_dir/$output.bam RGPL=$RGPL RGPU=$file_name RGSM=$file_name RGLB=$file_name 
    # samtools stuff
done < $1

#Run functions 
sort_bam
mark_duplicates
#consensus_sequence FIX_ME
# need to find a program to take fastq to fasta ( or write on myself ) should be easy

#coverage_plots
make_coverage_plots

call_variants
vcf_to_haplogrep

# rmdup is the filename
#if [[ $variant_caller == "gatk" ]]; then
#    echo "Gatk Variant Caller"
#    java \
#     -jar GenomeAnalysisTK.jar \
#     -T HaplotypeCaller \
#     -R $rCRS.fa \
#     `for i in ${tmp_dir}/*.rmdup.bam
#     do
#        echo "-I ${i}"
#     done`\
#     --validation_strictness LENIENT\
#     -o $tmp_dir/haplogrep.vcf

#elif [[ $variant_caller == "freebayes" ]];then
#    echo "Freebayes run :)"
#fi

