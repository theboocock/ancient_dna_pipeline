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
#PATH=$PATH:$DIR/../bin
#DIR=$( dirname "$(readlink  $0)")

#input setup file 

sort_bam(){    
parallel -j ${CORES} "java ${XMX} -jar ${PICARD}/SortSam.jar INPUT={} OUTPUT={.}.sorted.bam SORT_ORDER=coordinate" ::: $SAM_SEARCH_EXPAND
    #parallel -j 6 "java -Xmx1g -jar ${PICARD}/ValidateSamFile.jar INPUT={} OUTPUT=${tmp_dir}/{/.}.final.txt" ::: ${tmp_dir}/*.rmdup.bam
    SAM_SEARCH_EXPAND=${tmp_dir}/*.sorted.bam
}
#MOVE ME
mark_duplicates(){
    parallel -j ${CORES} "java ${XMX} -jar ${PICARD}/MarkDuplicates.jar INPUT={} OUTPUT=${tmp_dir}/{/.}.rmdup.bam METRICS_FILE=${tmp_dir}/{/.}.mark_dups.log AS=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT" ::: $SAM_SEARCH_EXPAND
}

SAM_SEARCH_EXPAND=${tmp_dir}/*.sorted.rmdup.bam

consensus_sequences(){
    parallel -j ${CORES} "samtools mpileup -uf ${reference} {} | bcftools view -cg - | vcfutils.pl vcf2fq | sed '/^@/!d;s//>/;N' > ${results_dir}/{/.}.fa" ::: $SAM_SEARCH_EXPAND
}

index_bams(){
    parallel -j ${CORES} "java ${XMX} -jar ${PICARD}/BuildBamIndex.jar \
        INPUT={}" ::: $SAM_SEARCH_EXPAND
}
call_variants_samtools(){
    vcf_output=$SETUP_FILE.vcf
    samtools mpileup -uf ${reference}.fa ${tmp_dir}/*.rmdup.bam  | bcftools view -cg - | vcfutils.pl varFilter -d ${MIN_DEPTH}> $results_dir/$vcf_output 
}

haplotype_caller(){
    #ploidy 100 for mtDNA
    # ss vcf is the single sample VCF
    parallel -j ${CORES} "${JAVA7} -jar ${XMX} ${GATK} \
        -T HaplotypeCaller \
        -R ${reference} \
        -I {} \
        -o ${tmp_dir}/{/.}.ss.vcf \
        -stand_emit_conf 10.0 " ::: $SAM_SEARCH_EXPAND
}
haplocaller_combine(){
    parallel -j ${CORES} "${JAVA7} -jar $XMX $GATK \
        -T HaplotypeCaller \
        -R ${reference} \
        --emitRefConfidence GVCF --variant_index_type LINEAR \
        --variant_index_parameter 128000 \
        -I {} \
        -o ${tmp_dir}/{/.}.gvcf" ::: $SAM_SEARCH_EXPAND
    $JAVA7 -jar -Xmx4g $GATK \
    -T GenotypeGVCFs \
    -R ${reference} \
    `for i in ${tmp_dir}/*.gvcf
    do
        echo "--variant ${i}"
    done` \
    -o ${vcf_output} 
}
vcf_output=$results_dir/$SETUP_FILE.raw.vcf

vcf_filter(){
   vcf_input=$vcf_output 
   $JAVA7 -jar $GATK \
    -T VariantFiltration \
    -R $reference \
    -V $vcf_output \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filterName "my_snp_filter" \
    -o $results_dir/$SETUP_FILE.filter.vcf 
}

generate_dna_damage_plots(){
  parallel -j ${CORES}  "mapDamage -i {} -d ${results_dir}/damage/{/.}" ::: ${SAM_SEARCH_EXPAND}
}

vcf_to_haplogrep(){
    # TODO add back the single HSD files although I feel they are not the usefuly
    #parallel -j ${CORES} "vcf_to_haplogrep.py -i {} -o ${results_dir}/{/.}.hsd" ::: ${tmp_dir}/*.ss.vcf
    vcf_to_haplogrep.py -i ${vcf_output} -o ${results_dir}/final_haplo.hsd
    vcf_to_haplogrep.py -i $results_dir/$SETUP_FILE.filter.vcf -o ${results_dir}/final_haplo_filtered.hsd
}

vcf_to_fasta(){
    parallel -j ${CORES} "vcf_to_fasta.py -i {} -o ${results_dir}/{/.}.fa -r ${reference} -s human" ::: ${tmp_dir}/*.ss.vcf
    vcf_to_fasta.py -i ${vcf_output} -o ${results_dir}/final_fasta_filtered.fa -r ${reference} 
}

#make_coverage_plots.py has hardcoded settings for coverage
coverage_plots(){
    parallel -j ${CORES} "samtools mpileup -D {} > $tmp_dir/{/.}.cov" ::: ${tmp_dir}/*.rmdup.bam
    make_coverage_plots.pl -- $tmp_dir/*.cov > $results_dir/coverage_plots.ps

}
coverage_plots_R(){
    parallel -j ${CORES} "samtools mpileup -D {} > $tmp_dir/{/.}.tmpcov" ::: ${tmp_dir}/*.rmdup.bam
    parallel -j ${CORES} "cat {} | cut -d $'\t' -f 1,2,3,4 > $tmp_dir/{/.}.cov" ::: ${tmp_dir}/*.tmpcov
    parallel -j ${CORES} "$RSCRIPTS/coverage_script.R -c {} -d ${results_dir}/coverage -s {/.} -o {/.} -r ${reference} -t ${results_dir}/coverage/coverage_data.txt" ::: $tmp_dir/*.cov 
}

add_and_or_replace_groups(){
    while read line
    do
        file_name=$(echo $line | cut -f -1 -d ' ')
        output=$file_name
        java ${XMX} -jar ${PICARD}/AddOrReplaceReadGroups.jar INPUT=${tmp_dir}/$output.sam OUTPUT=${tmp_dir}/$output.bam RGPL=$RGPL RGPU=$file_name RGSM=$file_name RGLB=$file_name 
done < $SETUP_FILE
SAM_SEARCH_EXPAND=${tmp_dir}/*.bam
}

map_reads(){
while read line
do
    file_name=$(echo $line | cut -f 1 -d ' ')
    output=$file_name 
    echo $reference
    if [[  $MAPPER =  bwa ]]; then
        if [[ $END = pe ]]; then
            # use bwa mem for paired end alignments
            pe_one=`echo ${line} | cut -d ' ' -f 2`
            pe_two=`echo ${line} | cut -d ' ' -f 3`
            bwa mem -t $CORES $reference ${pe_one} ${pe_two} > $tmp_dir/$output.sam 2> $tmp_dir/$output.bwa.err
        else
            se=$(echo $line | cut -f -d ' ' -f 2)
            bwa mem -t $CORES $reference $se > $tmp_dir/$ouput.sam 2> $tmp_dir/$ouput.bwa.err
        fi
            # use aln for se 
    elif [[ $MAPPER = "bowtie" ]]; then 
        # bowtie alignments
        if [[ $END = pe ]]; then
            bowtie2 -x $reference -p 6 -1 $pe_one -2 $pe_two > $tmp_dir/$output.sam 2> $tmp_dir/$output.bowtie2.err  
        else
            se=$(echo $line | cut -f -d ' ' -f 2)
            bowtie2 -x $reference -u $se > $tmp_dir/$output.sam 2> $tmp_dir/$output.bowtie.err
        fi 
    fi 
    echo $MAPPER
    # samtools stuff
done < $SETUP_FILE
SAM_SEARCH_EXPAND="${tmp_dir}/*.sam"
}
map_damage(){
while read line
do
    mapDamage
done < $SETUP_FILE
}

# $@ for add or replace groups
#Run functions 
#sort_bam
#mark_duplicates
#consensus_sequence FIX_ME
# recalibrated bam file
# need to find a program to take fastq to fasta ( or write on myself ) should be easy

#coverage_plots
#make_coverage_plots

# Don't use samtools
#call_variants_samtools
#vcf_to_haplogrep

# rmdup is the filename
#if [[ $variant_caller == "gatk" ]]; then
#    echo "Gatk Variant Caller"
#    java \
#     -jar GenomeAnalysisTK.jar \
#     -T HaplotypeCaller \
#     -R $reference.fa \
#     `for i in ${tmp_dir}/*.rmdup.bam
#     do
#        echo "-I ${i}"
#     done`\
#     --validation_strictness LENIENT\
#     -o $tmp_dir/haplogrep.vcf

#elif [[ $variant_caller == "freebayes" ]];then
#    echo "Freebayes run :)"
#fi

