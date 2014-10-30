import argparse

def read_length():

def main():
    parser = argparse.ArgumentParser(description="Filter by Read Length")
    parser.add_argument('-m',dest="min_length", help="Minimun_read_length")
    parser.add_argument('paired_end_reads',nargs="2",help="Fastq Input File")
    args = parser.parse_args()
    assert len(args.paired_end_reads) != 2, \
            "Require Paired end reads to be passed to file"
    assert args.min_length is not None, \
            "Require a minimum length for filtering"
    read_length(args.min_length,
if __name__=="__main__":
    main()
