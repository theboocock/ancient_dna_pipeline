#!/usr/bin/env python
#
# De interleave fastq file.
#

import argparse

def de_interleave(fastq)
    name_prefix = fastq.split('.fastq')[0]
    r1 = name_prefix + '_R1.fastq'
    r2 = name_prefix + '_R2.fastq'
    r1 = open(r1,'w')
    r2 = open(r2,'w')
    with open(fastq) as f:
        [r1.write(line) if (i % 8 < 4) else r2.write(line) for i, line in enumerate(open('test.fastq'))]
    r1.close()
    r2.close()

def main():
    parser=argparse.ArgumentParser(desciption="De-Interleavefastq")
    parser.add_argument('fastq',help="Interleaved Fastq")
    args = parser.parse_args()
    de_interleave(args.fastq)

if __name__=="__main__":
    main():
