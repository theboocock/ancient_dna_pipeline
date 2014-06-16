    echo "java \
     -Xmx4g -jar ~/bioinformatics/ancient_dna_pipeline/src/gatk/GenomeAnalysisTK.jar \
     -T HaplotypeCaller \
     -R $RCRS.fa \
     `for i in results/*.rmdup.bam
     do
        echo "-I ${i}"
     done`\
     --validation_strictness LENIENT\
     -o output.raw.snps.indels.vcf"
