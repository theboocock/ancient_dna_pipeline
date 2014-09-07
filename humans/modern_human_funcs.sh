# Author James Boocock
# Runs the pipeline on modern human populations
#
# 
# 2/14

indel_alignment(){
    # Indel alignment when I return
    INTERVALS=$tmp_dir/intervals.tmp
    ${JAVA7} ${XMX} ${GATK} \
    -T RealignerTargetCreator \
    -R ${reference} \
    -o tmp.intervals \
    -known $INDEL_REALIGN
    parallel -j ${CORES} "java $XMX -jar $GATK \ 
        -I {}
        -R ${reference}"
}

recalibrate_bam(){
    parallel -j ${CORES} "${JAVA7} -jar $XMX $GATK \
        -T PrintReads \
        -R ${reference} \
        -i {} \
        -BQSR {/.}.recal_rprt.grp \
        -o {/.}.recal.bam
        " ::: $SAM_SEARCH_EXPAND
        SAM_SEARCH_EXPAND=$tmp_dir/*.recal.bam
}

