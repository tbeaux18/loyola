#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 01/15/2019

find_motif.py

Input:
    Standard Input

"""

import sys


def generate_kmer(dna_sequence, kmer_length=3):
    """ Takes a DNA sequence and generates all kmer substrings of the sequence
        Args:
            dna_sequence (str) : must be ACGT; no RNA transcripts
            kmer_length (int) : default is set to 3 for codon, but can be changed
        Returns:
            a list of all possible kmers
    """

    kmer_list = [] # initializes the kmer list; order is important

    # iterates through the sequence creating a kmer length sliding window
    for i in range(len(dna_sequence) - kmer_length + 1):
        # appends each kmer to the kmer list
        kmer_list.append(dna_sequence[i:i+kmer_length])

    return kmer_list


def find_all_motifs(sequence, sub_sequence):
    """ Takes a sequence and its subsequence and finds all 1-based index find_all_motifs
        Args:
            sequence (str) : main dna sequence
            sub_sequence (str) : substring
        Returns:
            tab joined str of integers of all locations
    """

    list_of_kmers = generate_kmer(sequence, kmer_length=len(sub_sequence))

    result = [i + 1 for i, seq in enumerate(list_of_kmers) if seq == sub_sequence]

    return "\t".join(str(i) for i in result)


def main():
    """ takes std input, finds all motifs using 1-based indexing, writes to
        location_output.txt
    """

    std_input = [x.rstrip() for x in sys.stdin]

    sequence = std_input[0]
    sub_seq = std_input[1]

    with open('locations.txt', 'w') as location_output:
        location_output.write(find_all_motifs(sequence, sub_seq))


if __name__ == '__main__':
    main()
