#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 01/16/2019

shared_motif.py

Input:
    FASTA format txt/fa file


"""

import sys
from Bio import SeqIO
from suffix_trees import STree



def main():
    """ runs main """

    # requires the fasta/txt file is in the same directory
    fasta_text = sys.argv[1]

    fasta_list = list(SeqIO.parse(fasta_text, "fasta"))

    sequence_list = []

    for record in fasta_list:
        sequence_list.append(str(record.seq))

    st = STree.STree(sequence_list)

    with open('output-lcs.txt', 'w') as output:
        output.write(st.lcs())


if __name__ == '__main__':
    main()
