#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-18-2019

3A.py

A) (3 pts) Implement a dynamic programming algorithm for finding the length of a longest path in a
matrix given the following pseudocode. In this pseudocode, downi,j and righti,j are the respective weights
of the vertical and horizontal edges entering node (i,j). We denote the matrices holding (downi,j) and
(righti,j) as DOWN and RIGHT, respectively.

Input: (Omit variable names, just have integer inputs)
    m_length \n
    n_length \n
    right weight matrix equal to the value of m + 1 in an n-d array \n
    down weight matrix equal to the value of n + 1 in an n-d array \n

Output:
    max score of the last node

Example Test Case:
    # m_length = 4
    # n_length = 4
    # right = [[3, 2, 4, 0], [3, 2, 4, 2], [0, 7, 3, 4], [3, 3, 0, 2], [1, 3, 2, 2]]
    # down = [[1, 0, 2, 4, 3], [4, 6, 5, 2, 1], [4, 4, 5, 2, 1], [5, 6, 8, 5, 3]]

    result = 34
"""
import sys
from numpy import zeros

def longest_path_manhattan(m, n, down, right):
    """ a dynamic programmining implementation of the manhattan tour problem
        Args:
            m (int) : horizontal rows
            n (int) : vertical columns
            down (2d lst) : outer array must contain n + 1 arrays, with inner elements
                            of each array equal to n
            right (2d lst) : outer array must contain m + 1 arrays, with inner elements
                            of each array equal to m
        Returns:
            score[n][m] (int) : max score at the nth,mth coordinate of the matrix
    """
    # initializes numpy zero matrix of int
    s_matrix = zeros((n+1, m+1), dtype=int)

    # initialize S[0,0] = 0
    s_matrix[0][0] = 0

    # initializes each node vertically according to its weighted edge
    for i in range(1, n+1):
        s_matrix[i][0] = s_matrix[i-1][0] + down[i-1][0]

    # initializes each node horizontally according to its weighted edge
    for j in range(1, m+1):
        s_matrix[0][j] = s_matrix[0][j-1] + right[0][j-1]

    # begins the iteration, searching for the max score between the down
    # and right movements, then returns the summed max score at the last node
    for i in range(1, n+1):
        for j in range(1, m+1):
            s_matrix[i][j] = max(s_matrix[i-1][j] + down[i-1][j], s_matrix[i][j-1] + right[i][j-1])

    # returns the score at the last coordinate
    return s_matrix[n][m]





def main():
    """ runs main script """

    # std input must be exact for program to work
    std_input = sys.stdin.read().splitlines()

    m_length = int(std_input[0])
    n_length = int(std_input[1])

    # begins to do input handling. Problem is that 2d array is read in as a string
    # must do a lot of string manipulation to reconstruct sequence
    right = std_input[2].split('], [')
    down = std_input[3].split('], [')

    right_new = [r.replace('[[', '').replace(']]', '').split(', ') for r in right]
    right_int = [[int(r_int) for r_int in right_line] for right_line in right_new]

    down_new = [d.replace('[[', '').replace(']]', '').split(', ') for d in down]
    down_int = [[int(d_int) for d_int in down_line] for down_line in down_new]

    # runs and stores result as an integer
    path_score = longest_path_manhattan(m_length, n_length, down_int, right_int)

    print(path_score)

if __name__ == '__main__':
    main()
