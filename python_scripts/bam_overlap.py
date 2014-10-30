#!/usr/bin/env python
#
# Get the overlapping reads in the bam file and create a newbam record containing them
#
# @date 28/11/14
# @author James Boocock
#


import argparse


def get_fastq_line(fastq_file):
    fastq_line = {}
    fastq_line['header'] = fastq_file.readline.strip()
    fastq_line['seq'] = fastq_file.readline.strip()
    fastq_line['descrip'] = fastq_file.readline.strip()
    fastq_line['quality'] = fastq_file.readline.strip()
    return (fastq_line)

def average_quality(q1,q2):
    """ 
        Return the average of two quality scores
    """
    print str(unichr(((ord(q1) - 33) + (ord(q2) - 33))/ 2.0)) 

    return str(((ord(q1) - 33) + (ord(q2) - 33))/ 2.0)    

def rev_comp(item):
    """
        Retuns the base on the opposite strand
    """
    if item == 'A':
        return 'T'
    elif item == 'G':
        return 'C'
    elif item == 'T':
        return 'A'
    elif item == 'C':
        return 'T'
    else:
        return 'N'


def h(mat,r1seq,i,r2seq,j):
    """
    
    Function H for calculating the smith-waterman
    
    For now ignore gapped alignments.

    """
    print(r1seq)
    print(i)
    print(r2seq)
    print(j)
    match =mat[i-1][j-1] + s(r1seq[i-1],r2seq[j-1])
    return(max(match,0))    

def s(a,b):
    """ 
        Similarity matrix

    """
    if a == b:
        return 2
    elif a != b:
        return -1 

def smith_waterman(r1seq,r1qual,r2seq,r2qual):
    """ 
        Implement the smith waterman algorithm
        for sequence alignment.

        http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm

        Everyone should write one suppsoably.
        >>> smith_waterman('AT','TA')
        'AT'
    """
    m = len(r1seq)
    n = len(r2seq)
    #r2seq = reverse_complement(r2seq)
    mat = [[0 for x in xrange(n + 1)] for x in xrange(m+ 1)]
    maximum = [0,0,0] 
    for i in range(1,m+1):
        for j in range(1,n+1):
           mat[i][j] = h(mat, r1seq, i, r2seq, j)
           if ( mat[i][j] > maximum[0] ):
               maximum=[mat[i][j], i, j]
    i= maximum[1]
    j = maximum[2]
    r1_i = []
    r2_i = []
    rev_seq = '' 
    rev_qual = ''
    print(r1qual)
    print(r2qual)
    while(mat[i][j] != 0):
        rev_seq += r1seq[i-1]
        rev_qual+= average_quality(r1qual[i-1],r2qual[j-1]) 
        i = i - 1
        j = j-  1

    return(mat)

def reverse_complement(r2seq):
    """
        Returns the reverse complement of the second reads 
        IMPORTANT!
        >>> reverse_complement('CGGGG')
        'CCCCG'
        >>> reverse_complement('A')
        'T'
        >>> reverse_complement('ATATG')
        'CATAT'
    """
    return(("".join([rev_comp(item) for item in list(r2seq)]))[::-1])


def merge_reads(fastq_files):
    """ 
        Merge the reads and generate a consensus. 
    """
    with open(fastq_files[0]) as read_1:
        with open(fastq_files[1]) as read_2:
            r1 = get_fastq_line(read_1)
            r2 = get_fastq_line(read_2)
            r1seq = r1['seq']
            r2seq = r2['seq']
    
               
        
        

def main():
    parser = argparse.ArgumentParser(description="Get Overlapping Reads")
    parser.add_argument('fastq_files',nargs=2,help="Both pairs of the reads")
    args=parser.parse_args()
    merge_reads(fastq_files)
    


if __name__=="__main__":
    main()
