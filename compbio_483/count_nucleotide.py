#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 10/10/2018

Counting DNA Nucleotides Rosalind Problem 7
"""

import sys


def count_nucleotide(dna_sequence):
    """ Takes a sequence object and counts ACGT.
        Args:
            dna_sequence (str) : must contain ACGT
        Returns:
            ACGT integer counts in lexographical order
    """

    a_count, c_count, g_count, t_count = 0, 0, 0, 0 #initialize all counts to 0

    for nucleotide in dna_sequence.upper(): #iterate through str sequence
        #adds 1 to a count if nucleotide is true
        if nucleotide == 'A':
            a_count += 1
        elif nucleotide == 'C':
            c_count += 1
        elif nucleotide == 'G':
            g_count += 1
        elif nucleotide == 'T':
            t_count += 1

    # integer replacement to format output style
    return "%d %d %d %d" % (a_count, c_count, g_count, t_count)


def main():
    """ Takes sequence object and prints the counts in ACGT to console """
    seq_line = sys.stdin.read() #standard input
    print(count_nucleotide(seq_line))



if __name__ == '__main__':
    main()
