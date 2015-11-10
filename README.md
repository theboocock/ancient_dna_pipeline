ancient DNA pipeline.
====================

Processes ancient DNA data that was generated using NGS.

This pipeline employs a best practices approach for taking the data
from Raw FASTQ files and creates a number of useful output files.


Installation
============

Run ```./install.sh``` but before running make sure the INSTALL_DIR variable in this script points to a directory
that is located on your path.

#### Python

Python dependencies can be installed by running.

    pip install -r requirements.txt

This will install pysam, PyVCF, pyfasta, and Biopython.

#### mapDamage

Navigate to ```src/mapDamage/``` and follow all installation instructions at http://ginolhac.github.io/mapDamage/.

#### SeqMagick

Navigate to ```src/seqmagick/``` and run ```python setup.py install```

#### AdapterRemoval

Navigate to ```src/AdapterRemoval/``` and run the following commands.

    tar xvzf AdapterRemoval-1.5.4.tar.gz
    make

Then ensure that the executable ```AdapterRemoval``` is on your path. 
#### External dependencies.

All the following executable must be accesible from your path.

- Muscle (http://www.drive5.com/muscle/)
- R programming language (https://www.r-project.org/) with the packages
    - ggplot2, Biostrings, getopt installed
- bwa (http://bio-bwa.sourceforge.net/)
- samtools (http://samtools.github.io/)
- bcftools (https://samtools.github.io/bcftools/)


#### Example Run of the pipeline.

Navigate into the ```tests/test_data/``` directory and run the following commands.

    # generate pipeline_setup.txt 
    create_pipeline_setup.sh . raw
    # run pipeline for human mtDNA. make sure you replace the path to the reference file. 
    ancient_pipeline.sh -C "gi|251831106|ref|NC_012920.1|" -r ~/Programming/OpenSource/MyGitHub/ancient_dna_pipeline/ref/contamination.fa  \ -S "human" -t  -P 1
    







