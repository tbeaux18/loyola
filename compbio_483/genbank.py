#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-07-2019

genbank.py

"""

import sys
from Bio import Entrez


def main():
    """ runs main script """

    std_input = [x.rstrip() for x in sys.stdin]

    Entrez.email = "tbeaux18@gmail.com"

    handle = Entrez.esearch(db="nucleotide", \
    term=' "{}"[Organism] \
    AND ("{}"[PDAT] : "{}"[PDAT])'.format(std_input[0], std_input[1], std_input[2]))

    record = Entrez.read(handle)

    with open("genbank-output.txt", 'w') as output:
        output.write(str(record["Count"]))


if __name__ == '__main__':
    main()
