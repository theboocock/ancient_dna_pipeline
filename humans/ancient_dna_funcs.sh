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
annotate_traits(){
    parallel -j ${CORES} "make_traits.py -i {} -t $TRAITS_FILE -o ${results_dir}/{/.}.traits.nex" :::
}


sort_bam(){    
parallel -j ${CORES} "java ${XMX} -jar ${PICARD}/SortSam.jar INPUT={} OUTPUT={.}.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT" ::: $SAM_SEARCH_EXPAND
    #parallel -j 6 "java -Xmx1g -jar ${PICARD}/ValidateSamFile.jar INPUT={} OUTPUT=${tmp_dir}/{/.}.final.txt" ::: ${tmp_dir}/*.rmdup.bam
    SAM_SEARCH_EXPAND=${tmp_dir}/*.sorted.bam
}
#MOVE ME
mark_duplicates(){
    parallel -j ${CORES} "java ${XMX} -jar ${PICARD}/MarkDuplicates.jar INPUT={} OUTPUT=${tmp_dir}/{/.}.rmdup.bam METRICS_FILE=${tmp_dir}/{/.}.mark_dups.log AS=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT" ::: $SAM_SEARCH_EXPAND
SAM_SEARCH_EXPAND=${tmp_dir}/*.sorted.rmdup.bam
}


consensus_sequences(){
    parallel -j ${CORES} "samtools mpileup -uf ${reference} {} | bcftools view -cg - | vcfutils.pl vcf2fq | sed '/^@/!d;s//>/;N' > ${results_dir}/{/.}.fa" ::: $SAM_SEARCH_EXPAND
}

index_bams(){
    parallel -j ${CORES} "java ${XMX} -jar ${PICARD}/BuildBamIndex.jar \
        INPUT={} VALIDATION_STRINGENCY=LENIENT" ::: $SAM_SEARCH_EXPAND
}
call_variants_samtools(){
    vcf_output=$SETUP_FILE.vcf
    samtools mpileup -uf ${reference} ${tmp_dir}/*.rmdup.bam  | bcftools view -cg - | vcfutils.pl varFilter -d ${MIN_DEPTH}> $results_dir/$vcf_output 
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
    if [[ $BAIL_POS = "HAPLOCALLER" ]]; then
        echo "BROKE AT HAPLOTYPE_CALLER"
        exit 1
    fi
    vcf_output=$results_dir/$SETUP_FILE.raw.vcf
    if [[ $MAP_DAMAGE = "" ]]; then
        parallel -j ${CORES} "${JAVA7} -jar ${XMX} ${GATK} \
            -T HaplotypeCaller \
            -R ${reference} \
            --emitRefConfidence GVCF --variant_index_type LINEAR \
            --sample_ploidy $PLOIDY \
            --variant_index_parameter 128000 \
            -I {} \
            -o ${tmp_dir}/{/.}.gvcf" ::: $SAM_SEARCH_EXPAND
    elif [[ $MERGED_READS_ONLY = "TRUE" ]]; then
        echo "Using the second_read"
        ## TODO FIX this hard coding of paths.
        
        SAM_SEARCH_EXPAND="$results_dir/bams/*final_contaminant_mapped.bam"
        parallel -j ${CORES} "${JAVA7} -jar ${XMX} ${GATK} \
            -T HaplotypeCaller \
            -R ${reference} \
            --emitRefConfidence GVCF --variant_index_type LINEAR \
            --variant_index_parameter 128000 \
            --sample_ploidy $PLOIDY \
            -I {} \
            -o ${tmp_dir}/{/.}.gvcf" ::: $SAM_SEARCH_EXPAND
    else 
        echo "All the reads" 
        # This sholud hopefully fix this part of the analysis 
        # e.g running xapply and merging them both into the variant calling pipeline
        MERGED_READS_1=$results_dir/bams/*yes_collapse*.bam
        MERGED_READS_2=$results_dir/bams/*no_collapse*.bam 
        if [[ $TEST_FREEBAYES = "TRUE" ]]; then
            parallel -j ${CORES} --xapply "samtools merge {1} {2} > ${results_dir}/merged/{1/.}.merged.bam" ::: $MERGED_READS_1 ::: $MERGED_READS_2
            parallel -j ${CORES} "samtools index {}" ::: ${results_dir}/merged/*.bam
           freebayes -f ${reference} --ploidy 2 ${results_dir}/merged/*.bam > ${results_dir}/diploid_fb.vcf 
           freebayes -f ${reference} --ploidy 1 ${results_dir}/merged/*.bam > ${results_dir}/haploid_fb.vcf
        fi
        MERGED_READS_1=$results_dir/bams/*yes_collapse*.bam
       MERGED_READS_2=$results_dir/bams/*no_collapse*.bam 
       parallel -j 1 --xapply echo {1} {2} ::: $MERGED_READS_1 ::: $MERGED_READS_2
       parallel  -j ${CORES} --xapply "${JAVA7} -jar $XMX $GATK \
            -T HaplotypeCaller \
            -R ${reference} \
            --emitRefConfidence GVCF --variant_index_type LINEAR \
            --sample_ploidy $PLOIDY \
            --variant_index_parameter 128000 \
            -I {1} -I {2}  \
            --output_mode EMIT_ALL_SITES --allSitePLs \
            -o ${tmp_dir}/{1/.}.gvcf" ::: $MERGED_READS_1 ::: $MERGED_READS_2
    fi
    $JAVA7 -jar $XMX $GATK \
        -T GenotypeGVCFs \
        -R ${reference} \
        `for i in ${tmp_dir}/*.gvcf
do
    echo "--variant ${i}"
done` \
    --standard_min_confidence_threshold_for_calling 10 \
    -o ${vcf_output}
    $JAVA7 -jar $XMX $GATK \
        -T GenotypeGVCFs \
        -R ${reference} \
        `for i in ${tmp_dir}/*.gvcf
do
    echo "--variant ${i}"
done` \
    -o ${results_dir}/all_sites.vcf --includeNonVariantSites
}

vcf_filter(){
    vcf_input=$vcf_output 
    $JAVA7 $XMX -jar $GATK \
        -T VariantFiltration \
        -R $reference \
        -V $vcf_output \
        --filterExpression "MQ < 20.0 | DP < 3" \
        --filterName "my_snp_filter" \
        -o $results_dir/$SETUP_FILE.temp.vcf 
    cat $results_dir/$SETUP_FILE.temp.vcf | grep -v "my_snp_filter"  > $results_dir/$SETUP_FILE.filter.vcf
}
ancient_filter(){
    ancient_filter.py -d 2 --c2t --g2a -i {} -o $temp_results/{/.}.anc.fq
    
}
indel_realignment(){
    parallel --env PATH -j $CORES "$JAVA7 $XMX -jar $GATK \
        -T RealignerTargetCreator \
        -R $reference \
        -I {} \
        -o ${tmp_dir}/{/.}.intervals" ::: $SAM_SEARCH_EXPAND
    parallel --env PATH -j $CORES "$JAVA7 $XMX -jar $GATK \
        -T IndelRealigner \
        -R $reference \
        -I {}
        -targetIntervals ${tmp_dir}/{/.}.intervals
        -o {/.}.realigned.bam" ::: $SAM_SEARCH_EXPAND
        SAM_SEARCH_EXPAND=$tmp_dir/*.realigned.bam
}

map_damage(){
    parallel --env PATH -j ${CORES} "mapDamage -i {} -d ${results_dir}/damage/{/.} --rescale -r ${reference}" ::: ${SAM_SEARCH_EXPAND}
    parallel -j ${CORES} "cp {} ${tmp_dir}/" ::: ${results_dir}/damage/*/*bam
    SAM_SEARCH_EXPAND=${tmp_dir}/*rescaled.bam
}
pmd(){
    parallel --env PATH -j ${CORES} "samtools view -h {} | pmdtools.py --threshold 0 --header 2> {tmp_dir}/{/.}.pmd_stats.txt | samtools view -Sb - > ${tmp_dir}/{/.}.pmd_filter.bam"
    SAM_SEARCH_EXPAND=${tmp_dir}/*pmd_filter.bam
}

vcf_to_haplogrep(){
    # TODO add back the single HSD files although I feel they are not the usefuly
    #parallel -j ${CORES} "vcf_to_haplogrep.py -i {} -o ${results_dir}/{/.}.hsd" ::: ${tmp_dir}/*.ss.vcf
    if [[ $CONTAMINATION_MAPPING  = "" ]]; then
        vcf_output=$results_dir/$SETUP_FILE.raw.vcf
        vcf_to_fasta.py -i ${vcf_output} -o ${results_dir}/final_haplo.hsd -r ${reference} -s $SPECIES --use-indels --to-haplo --ploidy $PLOIDY
        create_table_haplogrep.py $vcf_output > ${results_dir}/final_haplo_table.txt
        vcf_to_fasta.py -i ${results_dir}/$SETUP_FILE.filter.vcf -o ${results_dir}/final_haplo_filtered.hsd -r ${reference} -s $SPECIES --use-indels --to-haplo --ploidy $PLOIDY  
        create_table_haplogrep.py $vcf_output > ${results_dir}/final_haplo_filtered_table.txt
    else
        vcf_output=$results_dir/$SETUP_FILE.raw.vcf
        vcf_to_fasta.py -i ${vcf_output} -o ${results_dir}/final_haplo.hsd -r ${reference} -s $SPECIES --use-indels --to-haplo --ploidy $PLOIDY -m $CONTAMINATION_MAPPING
        create_table_haplogrep.py $vcf_output > $results_dir/final_haplo_table.txt 
        vcf_to_fasta.py -i ${results_dir}/$SETUP_FILE.filter.vcf -o ${results_dir}/final_haplo_filtered.hsd -r ${reference} -s $SPECIES --use-indels --to-haplo --ploidy $PLOIDY -m $CONTAMINATION_MAPPING
        create_table_haplogrep.py ${vcf_output} > $results_dir/final_haplo_filtered.txt
    fi
}

vcf_to_snp_list(){
    # TODO add back the single HSD files although I feel they are not the usefuly
    #parallel -j ${CORES} "vcf_to_haplogrep.py -i {} -o ${results_dir}/{/.}.hsd" ::: ${tmp_dir}/*.ss.vcf
    if [[ $CONTAMINATION_MAPPING  = "" ]]; then
        vcf_output=$results_dir/$SETUP_FILE.raw.vcf
        vcf_to_fasta.py -i ${vcf_output} -o ${results_dir}/final_snp_list.txt -r ${reference} -s $SPECIES --use-indels --to-haplo --ploidy $PLOIDY --unique_only --to-haplo
        vcf_to_fasta.py -i ${results_dir}/$SETUP_FILE.filter.vcf -o ${results_dir}/final_snp_list_filtered.txt -r ${reference} -s $SPECIES --use-indels --to-haplo --ploidy $PLOIDY --unique_only --to-haplo
    if [[ $IMPUTATION = "TRUE" ]]; then
        vcf_to_fasta.py --ploidy $PLOIDY -i $impute_vcf -o $results_dir/final_snp_list_impute.txt -r $reference  -s $SPECIES --ploidy $PLOIDY $results_dir/coverage/files/*cov --unique --impute --to-haplo

        vcf_to_fasta.py --ploidy $PLOIDY -i $impute_vcf -o $results_dir/final_indels_snp_list_impute.txt -r $reference  --use-indels -s $SPECIES --ploidy $PLOIDY $results_dir/coverage/files/*cov --impute --unique --to-haplo
        vcf_to_fasta.py --ploidy $PLOIDY -i $impute_vcf -o $results_dir/final_indels_snp_list_impute_all.txt -r $reference  --use-indels -s $SPECIES --ploidy $PLOIDY $results_dir/coverage/files/*cov --impute  --to-haplo
    fi
    else
        vcf_output=$results_dir/$SETUP_FILE.raw.vcf
        vcf_to_fasta.py -i ${vcf_output} -o ${results_dir}/final_snp_list.txt -r ${reference} -s $SPECIES --use-indels --to-haplo --ploidy $PLOIDY -m $CONTAMINATION_MAPPING
        vcf_to_fasta.py -i ${results_dir}/$SETUP_FILE.filter.vcf -o ${results_dir}/final_snp_list_filtered.txt -r ${reference} -s $SPECIES --use-indels --to-haplo --ploidy $PLOIDY -m $CONTAMINATION_MAPPING

    if [[ $IMPUTATION = "TRUE" ]]; then
        vcf_to_fasta.py --ploidy $PLOIDY -i $impute_vcf -o $results_dir/final_snp_list_impute.txt -r $reference  -s $SPECIES --ploidy $PLOIDY $results_dir/coverage/files/*cov --unique --impute -m $CONTAMINATION_MAPPING --to-haplo

        vcf_to_fasta.py --ploidy $PLOIDY -i $impute_vcf -o $results_dir/final_indels_snp_list_impute.txt -r $reference  --use-indels -s $SPECIES --ploidy $PLOIDY $results_dir/coverage/files/*cov --impute --unique -m $CONTAMINATION_MAPPING --to-haplo
        vcf_to_fasta.py --ploidy $PLOIDY -i $impute_vcf -o $results_dir/final_indels_snp_list_impute_all.txt -r $reference  --use-indels -s $SPECIES --ploidy $PLOIDY $results_dir/coverage/files/*cov --impute  --to-haplo
    fi
fi
}
vcf_to_fasta(){
#    parallel -j ${CORES} "vcf_to_fasta.py -i {} -o ${results_dir}/{/.}.fa -r ${reference} -s human" ::: ${tmp_dir}/*.ss.vcf
    if [[ $CONTAMINATION_MAPPING != "" ]]; then
        vcf_to_fasta.py -i ${vcf_output} -o ${results_dir}/final_fasta.fa -r ${reference} -s $SPECIES --ploidy $PLOIDY   -m $CONTAMINATION_MAPPING $results_dir/coverage/files/*cov
        vcf_to_fasta.py -i ${results_dir}/$SETUP_FILE.filter.vcf -o ${results_dir}/final_fasta_filtered.fa -r ${reference} -s $SPECIES --ploidy $PLOIDY  -m $CONTAMINATION_MAPPING $results_dir/coverage/files/*cov

        # Add vcf_to_fasta
        vcf_to_fasta.py --ploidy $PLOIDY -i ${vcf_output} -o ${results_dir}/final_fasta_indels.fa -r ${reference} --use-indels -s $SPECIES --ploidy $PLOIDY  -m $CONTAMINATION_MAPPING $results_dir/coverage/files/*cov

        vcf_to_fasta.py --ploidy $PLOIDY -i ${results_dir}/$SETUP_FILE.filter.vcf -o ${results_dir}/final_fasta_filtered_indels.fa -r ${reference} --use-indels -s $SPECIES --ploidy $PLOIDY  -m $CONTAMINATION_MAPPING $results_dir/coverage/files/*cov

        if [[ $IMPUTATION = "TRUE" ]]; then
            vcf_to_fasta.py --ploidy $PLOIDY -i $impute_vcf -o $results_dir/final_fasta_impute.fa -r $reference  -s $SPECIES --ploidy $PLOIDY  -m $CONTAMINATION_MAPPING $results_dir/coverage/files/*cov --impute

            vcf_to_fasta.py --ploidy $PLOIDY -i $impute_vcf -o $results_dir/final_fasta_impute_indels.fa -r $reference  --use-indels -s $SPECIES --ploidy $PLOIDY  -m $CONTAMINATION_MAPPING $results_dir/coverage/files/*cov --impute
        fi
    else
        vcf_to_fasta.py -i ${vcf_output} -o ${results_dir}/final_fasta.fa -r ${reference} -s $SPECIES --ploidy $PLOIDY $results_dir/coverage/files/*cov
        vcf_to_fasta.py -i ${results_dir}/$SETUP_FILE.filter.vcf -o ${results_dir}/final_fasta_filtered.fa -r ${reference} -s $SPECIES --ploidy $PLOIDY $results_dir/coverage/files/*cov

        # Add vcf_to_fasta
        vcf_to_fasta.py --ploidy $PLOIDY -i ${vcf_output} -o ${results_dir}/final_fasta_indels.fa -r ${reference} --use-indels -s $SPECIES --ploidy $PLOIDY $results_dir/coverage/files/*cov

        vcf_to_fasta.py --ploidy $PLOIDY -i ${results_dir}/$SETUP_FILE.filter.vcf -o ${results_dir}/final_fasta_filtered_indels.fa -r ${reference} --use-indels -s $SPECIES --ploidy $PLOIDY $results_dir/coverage/files/*cov

        if [[ $IMPUTATION == "TRUE" ]]; then
            vcf_to_fasta.py --ploidy $PLOIDY -i $impute_vcf -o $results_dir/final_fasta_impute.fa -r $reference  -s $SPECIES --ploidy $PLOIDY $results_dir/coverage/files/*cov --impute

            vcf_to_fasta.py --ploidy $PLOIDY -i $impute_vcf -o $results_dir/final_fasta_impute_indels.fa -r $reference  --use-indels -s $SPECIES --ploidy $PLOIDY $results_dir/coverage/files/*cov --impute
        fi
    fi
}

fasta_to_nexus(){
    # Use python script to convert to nexus, removing all 
   seqmagick convert --output-format nexus --alphabet dna ${results_dir}/final_fasta.fa  $results_dir/temp.nex
   cat $results_dir/temp.nex | tr -d "'" > $results_dir/final.nex  
   seqmagick convert --output-format nexus --alphabet dna ${results_dir}/final_fasta_filtered.fa ${results_dir}/temp.nex 
   cat $results_dir/temp.nex | tr -d "'"  > $results_dir/final_filter.nex 
   seqmagick convert --output-format nexus --alphabet dna ${results_dir}/mus_fasta_filt.fa  $results_dir/temp.nex
   cat $results_dir/temp.nex | tr -d "'" > $results_dir/mus_filt.nex  
   seqmagick convert --output-format nexus --alphabet dna ${results_dir}/mus_fasta.fa ${results_dir}/temp.nex 
   cat $results_dir/temp.nex | tr -d "'"  > $results_dir/mus.nex 

}

align_muscle(){
    muscle -in $results_dir/final_fasta_filtered_indels.fa -out $results_dir/mus_fasta_filt.fa
    muscle -in $results_dir/final_fasta_indels.fa -out $results_dir/mus_fasta.fa
}
#make_coverage_plots.py has hardcoded settings for coverage

coverage_plots(){
    rm coverage/coverage_data.txt
    parallel -j ${CORES} "samtools mpileup -D {} > $tmp_dir/{/.}.cov" ::: ${tmp_dir}/*.rmdup.bam
    make_coverage_plots.pl -- $tmp_dir/*.cov > $results_dir/coverage_plots.ps

}

coverage_plots_R(){
    if [[ $MERGED_READS_ONLY = "" ]]; then
        MERGED_READS_1=$results_dir/bams/*yes_collapse*.bam
        MERGED_READS_2=$results_dir/bams/*no_collapse*.bam 
        parallel -j ${CORES} --xapply "samtools merge $results_dir/merged/{1/.}.merged.bam {1} {2} " ::: $MERGED_READS_1 ::: $MERGED_READS_2
        parallel -j ${CORES} "samtools index {}" ::: ${results_dir}/merged/*.bam
        SAM_SEARCH_EXPAND=${results_dir}/merged/*.bam
        parallel -j ${CORES} "samtools mpileup -D {} > $tmp_dir/{/.}.tmpcov" ::: ${results_dir}/merged/*.bam
    else
        parallel -j ${CORES} "samtools mpileup -D {} > $tmp_dir/{/.}.tmpcov" ::: ${SAM_SEARCH_EXPAND}
    fi
    parallel -j ${CORES} "cat {} | cut -d $'\t' -f 1,2,3,4 > $tmp_dir/{/.}.cov" ::: ${tmp_dir}/*.tmpcov
    mkdir ${results_dir}/coverage/files
    parallel -j ${CORES} " mv {} ${results_dir}/coverage/files/{/}" ::: $tmp_dir/*.cov
        parallel -j 1 "$RSCRIPTS/coverage_script.R -C \"${CONTAMINATION_MAPPING}\"   -d ${results_dir}/coverage -s {/.} -o {/.} -c {} -r ${reference} -t ${results_dir}/coverage/coverage_data.txt" ::: ${results_dir}/coverage/files/*.cov      
}

add_and_or_replace_groups(){
    while read line
    do
        file_name=$(echo $line | cut -f -1 -d ' ')
        output=$file_name
        if [[ $MAP_DAMAGE = "TRUE" ]]; then
            # Run mapDAmage on both the files
            # MS10148-2.no_collapse.rmdup.rescaled.sorted.sorted.bam 
            # MS10148-2.collapse.rmdup.rescaled.sorted.sorted.bam
            java ${XMX} -jar ${PICARD}/AddOrReplaceReadGroups.jar INPUT=${tmp_dir}/$output.collapse.rmdup.rescaled.ancient_filter.sorted.bam OUTPUT=${tmp_dir}/$output.yes_collapse.final.bam RGPL=$RGPL RGPU=$file_name RGSM=$file_name RGLB=$file_name  VALIDATION_STRINGENCY=LENIENT
        if [[ $MERGED_READS_ONLY = "" ]]; then
            java ${XMX} -jar ${PICARD}/AddOrReplaceReadGroups.jar INPUT=${tmp_dir}/$output.no_collapse.rmdup.rescaled.ancient_filter.sorted.bam OUTPUT=${tmp_dir}/$output.no_collapse.final.bam RGPL=$RGPL RGPU=$file_name RGSM=$file_name RGLB=$file_name  VALIDATION_STRINGENCY=LENIENT
        fi
        else
            java ${XMX} -jar ${PICARD}/AddOrReplaceReadGroups.jar INPUT=${tmp_dir}/$output.sorted.rmdup.bam OUTPUT=${tmp_dir}/$output.final.bam RGPL=$RGPL RGPU=$file_name RGSM=$file_name RGLB=$file_name VALIDATION_STRINGENCY=LENIENT 
        fi
    done < $SETUP_FILE
 SAM_SEARCH_EXPAND="${tmp_dir}/*.final.bam"
}

ancient_filter(){
   parallel -j $CORES "ancient_filter.py -i {} -o $temp_results/{/.}.ancient_filter.bam -f "bam" -g -c" :::  $SAM_SEARCH_EXPAND
   SAM_SEARCH_EXPAND=${tmp_dir}/*ancient_filter.bam
}

store_bams(){ 
    # TODO Problems to write to karen about.
    #if [[ $MERGED_READS_ONLY = ""]]; then 
    #    FIRST_READ=$(SAM_SEARCH_EXPAND | grep "yes_collapse")
     #   SECOND_READ=$(SAM_SEARCH_EXPAND | grep "no_collapse")

    #fi
    parallel -j $CORES "cp {} $results_dir/bams/{/.}.bam" ::: $SAM_SEARCH_EXPAND
    SAM_SEARCH_EXPAND="${results_dir}/bams/*.bam"
}
#recal_vcf(){
#    vcf_input=$results_dir/$SETUP_FILE.filter
#    if [[ $PLOIDY = "1" ]];
#        cat $results_dir/$SETUP_FILE..filter.vcf | sed "s/\t0:/\t0\/0:/g"  | sed "s/\t1:/\t1\/1:/g" | sed "s/\t2:/\t2\/2:/" | sed "s/\t\.:/\t\.\/\.:/g" > diploid.vcf
#    fi 
#}

beagle_imputation(){
    vcf_input=$results_dir/$SETUP_FILE.filter.recal.vcf  
    if [[ $PLOIDY = "1" ]]; then
        cat $results_dir/$SETUP_FILE.filter.vcf | perl -pe  "s/\t0:/\t0\/0:/g"  | perl -pe "s/\t1:/\t1\/1:/g" | perl -pe "s/\t2:/\t2\/2:/" | perl -pe "s/\t\.:/\t\.\/\.:/g" > $vcf_input 
        zero_info.py $vcf_input > tmp.tmp
        mv tmp.tmp $vcf_input
        percentage_imputed.py $vcf_input | sort -k 1 > ${results_dir}/imputed_percentage.txt
    elif [[ PLOIDY = "2" ]]; then
        recal_vcf.py $results_dir/$SETUP_FILE.filter.vcf >  $vcf_input
    fi 
    # This function crates the imputed files for further processing.
    java -jar $BEAGLE gt=${vcf_input} out=tmp
    gzcat tmp.vcf.gz > $results_dir/$SETUP_FILE.impute.vcf
    impute_vcf=$results_dir/$SETUP_FILE.impute.vcf
}

#remove_g_a_c_t(){
   #cat $results_dir/pipeline_setup.txt.filter.vcf | ack -v  "G\tA|C\tT" > $results_dir/ancient_strict.vcf
   # vcf_to_fasta.py -i $results_dir/ancient_strict.vcf -o ${results_dir}/ancient_strict.fa -r ${reference}
   #cat $results_dir/temp.nex | tr -d "'" > $results_dir/ancient_strict.nex
   #rm $results_dir/temp.nex
   #seqmagick convert --output-format nexus --alphabet dna ${results_dir}/ancient_strict.fa  $results_dir/temp.nex
   
#}
map_reads(){
while read line
do
    file_name=$(echo $line | cut -f 1 -d ' ')
    output=$file_name
    extra_arg_bwa=""
    if [[ $MAP_DAMAGE != "" ]]; then
        extra_arg_bwa="-n 0.03 -o 2 -l 1024"
    fi 
    if [[  $MAPPER =  bwa ]]; then
        if [[ $END = pe ]]; then
            # use bwa mem for paired end alignment
            pe_one=$(echo ${line} | cut -d ' ' -f 2)
            pe_two=$(echo ${line} | cut -d ' ' -f 3)
            if [[ $MAP_DAMAGE != "" ]]; then
                echo $pe_one
                echo $pe_two
                AdapterRemoval --collapse --file1 ${pe_one} --file2 ${pe_two} --outputstats ${pe_one}.stats --trimns --outputcollapsed $tmp_dir/${pe_one}.collapsed --minlength 25 --output1 $tmp_dir/$pe_one.p1 --output2 $tmp_dir/${pe_two}.p2 --mm 3  --minquality 20   --trimqualities
            fi
            if [[ $MAP_DAMAGE = "" ]]; then
            echo $pe_one
            echo $pe_two
            # -M is essensetion for picard compatibility
            bwa mem -M -t $CORES $reference ${pe_one} ${pe_two} > $tmp_dir/$output.sam 2> $tmp_dir/$output.bwa.err
            else
                bwa aln -t $CORES $extra_arg_bwa $reference $tmp_dir/$pe_one.collapsed > tmp1.sai
                bwa samse $reference tmp1.sai $tmp_dir/$pe_one.collapsed > $tmp_dir/$output.collapse.sam 2> $tmp_dir/$output.collapse.bwa.err
                if [[ $ONLY_MERGE_READS = "" ]]; then
                    bwa aln -t $CORES $extra_arg_bwa $reference $tmp_dir/$pe_one.p1 > tmp1.sai 
                    bwa aln -t $CORES $extra_arg_bwa $reference  $tmp_dir/$pe_two.p2 > tmp2.sai 
                    bwa sampe $reference tmp1.sai tmp2.sai $tmp_dir/$pe_one.p1 $tmp_dir/$pe_two.p2 > $tmp_dir/$output.no_collapse.sam 
                    samtools view -Sb $tmp_dir/$output.no_collapse.sam > $tmp_dir/$output.no_collapse.bam
                    samtools sort $tmp_dir/$output.no_collapse.bam $tmp_dir/$output.no_collapse.sorted
                # Mark Duplicates Paired END
                    java ${XMX} -jar ${PICARD}/MarkDuplicates.jar INPUT=${tmp_dir}/$output.no_collapse.sorted.bam OUTPUT=${tmp_dir}/$output.no_collapse.rmdup.bam METRICS_FILE=${tmp_dir}/$output.no_callapse.mark_dups.log AS=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT
                fi
                #FilterUniqueBam
                samtools view -Sb $tmp_dir/$output.collapse.sam > $tmp_dir/$output.collapse.bam
                samtools sort $tmp_dir/$output.collapse.bam $tmp_dir/$output.collapse.sorted
                filterUniqueBam.py --remove-duplicates < $tmp_dir/$output.collapse.sorted.bam > $tmp_dir/$output.collapse.rmdup.bam
            fi
        else
            se=$(echo $line | cut  -d ' ' -f 2)
            echo $se
            bwa mem -t $CORES ${extra_arg_bwa} $reference $se > $tmp_dir/$output.sam 2> $tmp_dir/$output.bwa.err
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
    # samtools stuff
done < $SETUP_FILE
SAM_SEARCH_EXPAND="${tmp_dir}/*.sam"

if [[ $MAP_DAMAGE != "" ]]; then
    SAM_SEARCH_EXPAND="${tmp_dir}/*.collapse.rmdup.bam"
    map_damage
    if [[ $MERGED_READS_ONLY = "" ]]; then
        SAM_SEARCH_EXPAND="${tmp_dir}/*.no_collapse.rmdup.bam"
        map_damage
    fi
    #merge these adjusted bam files
    # Downweight the reads in the final file. 
    # First 2 T's  = 2 quality
    # First 2 G's =  2 quality
    SAM_SEARCH_EXPAND="${tmp_dir}/*.rescaled.bam"
    parallel -j $CORES "ancient_filter.py -i {} -f \"bam\" -g -c -d 2 -o ${tmp_dir}/{/.}.ancient_filter.bam" ::: $SAM_SEARCH_EXPAND
    SAM_SEARCH_EXPAND="${tmp_dir}/*.rescaled.ancient_filter.bam"
fi
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

remove_contaminants(){
    parallel -j ${CORES} "samtools view -bh {} \"${CONTAMINATION_MAPPING}\" > ${tmp_dir}/{/.}_contaminant_mapped.bam" ::: $SAM_SEARCH_EXPAND
    SAM_SEARCH_EXPAND="${tmp_dir}/*_contaminant_mapped.bam"
}

save_contaminants(){
    parallel -j $CORES "cp {} $results_dir/contamination/{/.}.bam" ::: $SAM_SEARCH_EXPAND
    parallel -j $CORES "samtools index {} " ::: ${results_dir}/contamination/*.bam
}
contamination_percentage(){
    # Contamination,
    mkdir -p ${results_dir}/contamination_merged
    if [[ $MERGED_READS_ONLY = "" ]]; then
        MERGED_READS_1=$results_dir/contamination/*yes_collapse*.bam
        parallel -j ${CORES} --xapply "samtools merge $results_dir/contamination_merged/{1/.}.merged.bam {1} {2}" ::: $MERGED_READS_1 ::: $MERGED_READS_2
        parallel -j ${CORES} "samtools index {}" ::: ${results_dir}/contamination_merged/*.bam
        parallel -j $CORES "contamination_percentage.py -r \"${CONTAMINATION_MAPPING}\" {} > $results_dir/contamination_merged/{/.}.txt" :::  ${results_dir}/contamination_merged/*.bam 
    else
        MERGED_READS_1=$results_dir/contamination/*yes_collapse*.bam
        parallel -j $CORES "contamination_percentage.py -r \"${CONTAMINATION_MAPPING}\" {} > $results_dir/contamination_merged/{/.}.txt" :::  $MERGED_READS_1
    fi
   contamination_merged.py $results_dir/contamination_merged/*.txt > $results_dir/contamination_table.txt  
}

