ancient DNA pipeline.
====================

Processes ancient DNA data that was generated using NGS.

This pipeline employs a best practices approach for taking the data
from Raw FASTQ files and creates a number of useful output files.


Installation
============

Run ```./install.sh``` make sure the INSTALL_DIR variable in this script points to a directory
that is located on your path.

### Python

Python dependencies can be installed by running.

    pip install -r requirements.txt

This will install pysam, PyVCF, pyfasta, and Biopython.


### mapDamage


Navigate to ```src/mapDamage/``` and follow installation instructions at http://ginolhac.github.io/mapDamage/.

### SeqMagick

Navigate to ```src/seqmagick/``` and run ```python setup.py install```

### AdapterRemoval

Navigate to ```src/AdapterRemoval/``` and run the following commands.

    tar xvzf AdapterRemoval-1.5.4.tar.gz
    make

Then ensure that the executable ```AdapterRemoval``` is on your path. 
### External dependencies.

All the following executable must be accesible from your path.

- Muscle (http://www.drive5.com/muscle/)
- R programming language (https://www.r-project.org/) with the packages
    - ggplot2, Biostrings, getopt installed
- bwa (http://bio-bwa.sourceforge.net/)
- samtools (http://samtools.github.io/)
- bcftools (https://samtools.github.io/bcftools/)







