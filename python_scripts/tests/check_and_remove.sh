#!/usr/bin/bash

# Test check_and_remove
check_and_remove_samples.py -r $RCRS.fa test_files/MS10062.no_collapse.final_contaminant_mapped.bam -s "dog" -c "gi|17737322|ref|NC_002008.4|" -l 16727 -d 12
check_and_remove_samples.py -r $RCRS.fa test_files/MS10062.*.bam -s "dog" -c "gi|17737322|ref|NC_002008.4|" -l 16727 -d 12 -m

check_and_remove_samples.py -r $RCRS.fa test_files/MS10062.no_collapse.final_contaminant_mapped.bam -s "dog" -c "gi|17737322|ref|NC_002008.4|" -l 16727 -d 100
check_and_remove_samples.py -r $RCRS.fa test_files/MS10062.*.bam -s "dog" -c "gi|17737322|ref|NC_002008.4|" -l 16727 -d 100 -m

