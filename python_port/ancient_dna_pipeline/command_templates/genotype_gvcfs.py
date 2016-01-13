GENOTYPEGVCFS_TEMPLATE = """
    java -Djava.io.tmpdir=tmpdir -jar {0} \
    {1} \
    -T GenotypeGVCFs \
    -R {2} \
    -o {3} \
"""

