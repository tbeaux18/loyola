#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-11-2019

fastq2fasta.py

biopython as dependency

"""
import sys
from Bio import SeqIO


def main():
    """ runs main script """

    # not standard input
    # takes file name directly as argument
    fasta_text = sys.argv[1]

    # biopython to convert fastq to fasta
    # output name is fasta-output.txt
    SeqIO.convert(fasta_text, "fastq", 'fasta-output.txt', "fasta")


if __name__ == '__main__':
    main()
