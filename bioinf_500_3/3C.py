#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-19-2019

3C.py

Input:
    Sequence 1 \n
    Sequence 2 \n

"""

import sys
from numpy import zeros


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


def output_lcs(backtrack_matrix, v, i, j, lc_sequence):
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
            Horizontal sequence displays insertions

    """
    if i == 0 or j == 0:
        return

    if backtrack_matrix[i][j] == -1:
        output_lcs(backtrack_matrix, v, i - 1, j, lc_sequence)

    elif backtrack_matrix[i][j] == 2:
        output_lcs(backtrack_matrix, v, i, j-1, lc_sequence)
        lc_sequence.append('-')

    elif backtrack_matrix[i][j] == 1:
        output_lcs(backtrack_matrix, v, i - 1, j - 1, lc_sequence)
        lc_sequence.append(v[i - 1])


def main():
    """ runs main script """

    std_input = [x.strip() for x in sys.stdin]

    v_seq = std_input[0]
    w_seq = std_input[1]

    longest_common_subsequence = []

    backtrack_matrix = lcs_backtrack(v_seq, w_seq)

    output_lcs(backtrack_matrix, v_seq, len(v_seq), len(w_seq), longest_common_subsequence)

    joined_lcs = ''.join(longest_common_subsequence)

    print(joined_lcs)

if __name__ == '__main__':
    main()
