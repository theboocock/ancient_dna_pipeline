#!/usr/bin/env python
#
#
# First thing here is to create a spreadsheet with gaps as tabs for the sample that have missing information.
#
#
#
import itertools
import argparse
import re
from collections import OrderedDict

def tryint(s):
    try:
        return int(s)
    except:
        return s


def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def create_site_differences(haplogrep_file):
    with open(haplogrep_file) as h:
        snps = set()
        samples = OrderedDict()
        for line in h:
            line = line.strip().split('\t')
            sample = line[0]
            subs = line[3:]
            for sub in subs:
                snps.add(sub)
                try:
                    samples[sample].append(sub)
                except KeyError:
                    samples[sample] = [sub]
    #l = [["snp"],samples.keys()]
    #l = list(itertools.chain(*l))
    final_lines = [samples.keys()]
    for snp in sorted(snps,key=alphanum_key):
        out_line = [] 
        for sample, subs in samples.items():
            if snp in subs:
                out_line.append(snp)
            else:
                out_line.append(' ')
        final_lines.append(out_line)  
    for rows in zip(*final_lines):
        print '\t'.join(rows)


def main():
    parser = argparse.ArgumentParser(description="Takes a haplogrep file and outputs a nicely formatted tab-delimited file.")
    parser.add_argument('haplogrep_file',help='haplogrep_file')
    args = parser.parse_args()
    create_site_differences(args.haplogrep_file)
 
if __name__=="__main__":
    main()
