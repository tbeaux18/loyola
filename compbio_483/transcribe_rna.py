#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 01/15/2019

transcribe_rna.py

Input:
    standard Input of dna sequence
Output:
    creates output file rna_output.txt
"""

import sys


def main():
    """ transcribes dna to rna and prints to console """

    dna_input = "".join([x.rstrip() for x in sys.stdin])

    rna_output = dna_input.replace('T', 'U')

    with open('rna_output.txt', 'w') as output:
        output.write(rna_output)


if __name__ == '__main__':
    main()
