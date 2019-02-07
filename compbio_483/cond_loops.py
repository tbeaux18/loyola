#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 10/10/2018

Condition and Loops Rosalind Problem 4

"""

import sys


def count_odd_numbers(first_number, second_number):
    """
    Args:
        first_number (int) : starting number
        second_number (int) : end number
    Returns:
        Sum of all odd integers between first and second number
    """
    sum_total = 0 # initialize count to 0

    for number in range(first_number, second_number):
        if number % 2 == 1: # Checks divisibility to see if number is odd
            sum_total += number #adds the number to the sum_total

    return sum_total


def main():
    """ Takes 2 integers from stdin and prints total odd number count
    between the two integers
    """

    for line in sys.stdin: #iterates through stdin txt file
        line = line.split() #converts the first line from file into a list
        start_number = int(line[0]) #start number
        end_number = int(line[1]) #end number

    print(count_odd_numbers(start_number, end_number))



if __name__ == '__main__':
    main()
