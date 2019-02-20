#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-19-2019

3D.py

Input:
    Sequence 1 \n
    Sequence 2 \n
"""

import sys
from numpy import zeros


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



def lcs_backtrack(v, w):
    """ creates a trace back of two sequences of varying lengths

        Args:
            v (sequence) : can be any sequence, but object class must be same as
                        second sequence
            w (sequence) : can be any sequence, but object class must be same as
                        first sequence
        Returns:
            backtrack_array (2d array) : matrix where if:
                down = -1
                right = 2
                diagonal = 1

    """

    score = zeros((len(v)+1, len(w)+1), dtype=int)

    backtrack_array = zeros((len(v)+1, len(w)+1), dtype=int)

    for i in range(1, len(v) + 1):
        score[i][0] = 0

    for j in range(1, len(w) + 1):
        score[0][j] = 0

    for i in range(1, len(v) + 1):

        for j in range(1, len(w) + 1):

            if v[i - 1] == w[j - 1]:

                diag = score[i-1][j-1] + 1

            else:
                diag = score[i-1][j-1]

            score[i][j] = max(score[i-1][j], score[i][j-1], diag)

            if score[i][j] == score[i-1][j]: # down
                backtrack_array[i][j] = -1

            elif score[i][j] == score[i][j-1]: # right
                backtrack_array[i][j] = 2

            else: # score[i][j] == score[i-1][j-1] + 1 and v[i] == w[j]: # diagonal
                backtrack_array[i][j] = 1

    return backtrack_array


def build_output_lcs(backtrack_matrix, v, i, j, lc_sequence):
    """ concatenates the individual elements back to the longest
        common subsequence; recursively builds the string

        Args:
            backtrack_matrix (2d array) : backtrack matrix from lcs_backtrack
            v (str) : horizontal sequence
            i (int) : length of v
            j (int) : length of j
            lc_sequence (lst) : initialized empty list to build sequence
        Returns:
            Empty (not None) array is built outside the function
            Horizontal sequence displays edit distance
            Depending on length of the first or second string counts are made
            based on deletions or insertions
    """
    if i == 0 or j == 0:
        return

    #
    if backtrack_matrix[i][j] == -1:
        build_output_lcs(backtrack_matrix, v, i - 1, j, lc_sequence)
        lc_sequence.append('+')

    elif backtrack_matrix[i][j] == 2:
        build_output_lcs(backtrack_matrix, v, i, j-1, lc_sequence)
        lc_sequence.append('-')

    elif backtrack_matrix[i][j] == 1:
        build_output_lcs(backtrack_matrix, v, i - 1, j - 1, lc_sequence)
        lc_sequence.append(v[i - 1])


def main():
    """ runs main script """

    std_input = [x.strip() for x in sys.stdin]

    the_first_sequence = std_input[0]
    the_second_sequence = std_input[1]

    longest_common_subsequence = []

    backtrack_matrix = lcs_backtrack(the_first_sequence, the_second_sequence)

    build_output_lcs(backtrack_matrix, \
                the_first_sequence, \
                len(the_first_sequence), \
                len(the_second_sequence), \
                longest_common_subsequence)

    joined_lcs = ''.join(longest_common_subsequence)

    # the edit distance between 2 equal strings is the hamming distance
    if len(the_first_sequence) == len(the_second_sequence):
        print(hamming_distance(the_first_sequence, the_second_sequence))

    # if the strings are of unequal length, then depending on which string
    # was larger will determine which count of insertion or deletion to
    # perform on the overall sequence
    else:
        if len(the_first_sequence) > len(the_second_sequence):
            print(joined_lcs.count('+'))
        else:
            print(joined_lcs.count('-'))

if __name__ == '__main__':
    main()
