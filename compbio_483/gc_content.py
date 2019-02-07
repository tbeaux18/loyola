#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 01/15/2019

gc_content.py

Requirements/Dependencies:
    input file must be in FASTA format
    biopython for FASTA parser
Input: fasta txt file; must be in same directory; no external path accounted for.
Output: Title
        Absolute value of the sequence most GC content
"""

import sys
from Bio import SeqIO


def max_gc_count(fasta_seqio_list):
    """ takes in a list of SeqIO objects and returns the maximum
    value
        Args:
            fasta_seqio_list (list) : SeqIO objects
        Returns:
            >title + \n + max result
    """
    fasta_dict = {}

    for record in fasta_seqio_list:
        g_count = record.seq.count('G')
        c_count = record.seq.count('C')
        total_gc = (g_count + c_count) / len(record.seq) * 100
        fasta_dict[record.id] = total_gc

    max_result = max(fasta_dict, key=fasta_dict.get)

    return max_result + '\n' + str(fasta_dict[max_result])


def main():
    """ creates SeqIO objects, finds max gc count, and prints to console """

    # requires input file to be in same directory
    fasta_text = sys.argv[1]

    fasta_list = list(SeqIO.parse(fasta_text, "fasta"))

    print(max_gc_count(fasta_list))

if __name__ == '__main__':
    main()
