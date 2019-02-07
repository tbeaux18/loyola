#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 01/15/2019

trans_ratio.py

Requirements/Dependencies:
    Biopython
    txt file in fasta format
Input:
    FASTA format txt file in the same directory as this file; no external path
    accounted for.
Output:
    Ratio of transitions to transversions

"""

import sys
from Bio import SeqIO


def calculate_trans_ratio(sequence_one, sequence_two):
    """ Calculates transition/transversion ratio between 2 strings
        Args:
            sequence_one (str) : sequence of ACGT
            sequence_two (str) : sequence of ACGT
        Returns:
            trans_ratio (int) : number of mismatches between 2 strings
    """

    transition_count = 0
    transversion_count = 0

    for nucleotide_1, nucleotide_2 in zip(sequence_one, sequence_two):
        if nucleotide_1 != nucleotide_2:
            if nucleotide_1+nucleotide_2 in ('GA', 'AG', 'CT', 'TC'):
                transition_count += 1
            else:
                transversion_count += 1

    return transition_count / transversion_count




def main():
    """ creates SeqIO objects, stores sequences, runs the trans_ratio, prints to console"""

    # requires the fasta/txt file is in the same directory
    fasta_text = sys.argv[1]

    fasta_list = list(SeqIO.parse(fasta_text, "fasta"))

    sequence_tmp_list = []

    for record in fasta_list:
        sequence_tmp_list.append(record.seq)

    print(calculate_trans_ratio(sequence_tmp_list[0], sequence_tmp_list[1]))

if __name__ == '__main__':
    main()
