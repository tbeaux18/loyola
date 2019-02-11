#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-11-2019

2A.py

A) Generate the k-mer composition of a string.
Input: A string Text and an integer k.
Output: COMPOSITIONk(Text), where the overlapping k-mers in Text are written in lexicographic order.

Standard Input:
    kmer length (int) \n
    sequence (str) \n
        - NOTE: script will concatenate sequences if newline breaks are present
Output to txt file:
    lexographic kmer \n
"""


import sys

def composition_k(sequence, kmer_length):
    """ generates list of kmers of length k and sorts them
        lexographically
        Args:
            sequence (str)
            kmer_lenth (int)
        Returns:
            kmer_list (lst) : sorted lexographically
    """
    kmer_list = []

    for i in range(len(sequence) - kmer_length + 1):
        kmer_list.append(sequence[i:i+kmer_length])

    return sorted(kmer_list)



def main():
    """ takes std input and runs main script """

    std_input = sys.stdin.read().splitlines()

    kmer_lenth = int(std_input[0])
    comp_seq = std_input[1]

    with open('2A-output.txt', 'w') as output:
        # checks length of std_input, will concatenate subsequence strings
        # if length > 2 since newlines are assumed to be present
        if len(std_input) > 2:
            comp_seq = ''.join(std_input[1:])
            for kmer in composition_k(comp_seq, kmer_lenth):
                output.write("{}\n".format(kmer))
        else:
            for kmer in composition_k(comp_seq, kmer_lenth):
                output.write("{}\n".format(kmer))


if __name__ == '__main__':
    main()
