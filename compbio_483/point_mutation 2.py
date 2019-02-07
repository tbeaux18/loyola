#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 10/15/2018

Counting Point Mutations Rosalind Problem 10

"""

import sys


def hamming_distance(sequence_one, sequence_two):
    """ Calculates hamming distance between 2 strings
        Args:
            sequence_one (str) : sequence of ACGT
            sequence_two (str) : sequence of ACGT
        Returns:
            Hamming distance (int) : number of mismatches between 2 strings
    """

    count = 0 #initialize the count to 0
    for n_1, n_2 in zip(sequence_one, sequence_two): #iterates through 2 sequences
        if n_1 != n_2: #checks equality against each nucleotide
            count += 1 # if not equal adds 1 to the count
    return count



def main():
    """ Takes standard input of two sequences and prints the hamming hamming_distance
    to the console
    """

    #lines = [x.rstrip() for x in sys.stdin.readlines()] # converts std input into list
    first_sequence = 'AACCTTGGAATTT'
    second_sequence = 'ACCCTTGAGCTTT'


    print(hamming_distance(first_sequence, second_sequence))


if __name__ == '__main__':
    main()
