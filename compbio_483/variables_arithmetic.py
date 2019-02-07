#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 10/10/2018

Variables and Some Arithmetic Rosalind Problem 2

Input: Standard input
Output: prints function output to console

"""

import sys

def calculate_hyp(side_a, side_b):
    """ Returns the sum of the square of side a and side b.
    Args:
        side_a (int) : Side length of a right triangle
        side_b (int) : Side length of a right triangle

    Returns:
        hypotenuse_squared (int) : Sum of A^2 + B^2
    """
    # sums the square of side a and side b
    hypotenuse_squared = ((int(side_a) ** 2) + (int(side_b) ** 2))

    return int(hypotenuse_squared)


def main():
    """ prints function output to console """
    # reads in standard input file Ex. 3\n5
    for line in sys.stdin: #iterates through stdin txt file
        line = line.split() #converts the first line from file into a list
        length_a = line[0] #first position in list
        length_b = line[1] #second position in list

    print(calculate_hyp(length_a, length_b))



if __name__ == '__main__':
    main()
