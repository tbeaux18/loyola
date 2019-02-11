#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 10/10/2018

Complementing a Strand of DNA Rosalind Problem 8
"""

import sys

def generate_reverse_complement(dna_sequence):
    """
        Args:
            dna_sequence (str) : must be a string of ACGT
        Returns:
            Reverse complement of an ACGT string
    """

    complement = []

    for nucleotide in dna_sequence.upper(): #iterates through str dna_sequence

        # if conditions for complements and appends to complement list
        if nucleotide == 'A':
            complement.append('T')
        elif nucleotide == 'C':
            complement.append('G')
        elif nucleotide == 'G':
            complement.append('C')
        elif nucleotide == 'T':
            complement.append('A')


    return "".join(complement[::-1]) #joins list and returns reverse sequence




def main():
    """ Takes stdin and prints reverse comp to console """

    seq_line = sys.stdin.read().rstrip()

    print(generate_reverse_complement(seq_line))



if __name__ == '__main__':
    main()
