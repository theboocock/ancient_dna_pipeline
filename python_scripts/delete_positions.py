#!/usr/bin/env python
#
# @author James
# Remove  positions where C->T transitions occur
# in both the files.
#
# Only one use case for me currently, to delete SNPS that contain these tranisitions in the ancient Data
#
# So takes a muscle aligned file, then removes any positions that I have said to remove from the reference

import argparse



def main():
    parser = argparse.ArgumentParser(description="Muscle alignment remove ancient DNA SNPS from the final file")
    parser.add_argument('-i','--input-fa',dest='input_fa',help='Muscle aligned fasta input')
    parser.add_argument('-o','--offset',dest='offset', help='Starting position of the Ref Sequence')
    parser.add_argument('-r','--reference_name',dest='ref_seq', help="Sequence name of the reference sequence")
    parser.add_argument('-p','--positions',dest='position_file', help="Positions that could be aDNA damage")
    args = parser.parse_args()
     

if __name__=="__main__":
    main():
