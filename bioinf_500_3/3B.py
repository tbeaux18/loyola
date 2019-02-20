#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-19-2019

3B.py

Input:
    sequence 1 \n
    sequence 2 \n
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



def main():
    """ runs main script """

    std_input = [x.strip() for x in sys.stdin]

    v_seq = std_input[0]
    w_seq = std_input[1]

    backtrack_matrix = lcs_backtrack(v_seq, w_seq)

    print(backtrack_matrix)


if __name__ == '__main__':
    main()
