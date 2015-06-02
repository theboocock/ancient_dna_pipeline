#!/usr/bin/env python
#
# Merge bams  using pysam
#
# @author James Boocock
#


import argparse


def main():
    """
        Merges reads using the readgroup 
    """
    parser = argparse.ArgumentParser(description='Merge reads')
    args = parser.parse_args()


if __name__ == "__main__":
    main()
