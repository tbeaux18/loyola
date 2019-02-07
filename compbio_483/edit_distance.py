#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 01-30-2019

edit_distance.py

"""

import sys
import editdistance
from Bio import SeqIO


def main():
    """ runs main script """

    # requires the fasta/txt file is in the same directory
    fasta_text = sys.argv[1]

    fasta_list = list(SeqIO.parse(fasta_text, "fasta"))

    fasta_sequence = []

    for record in fasta_list:
        fasta_sequence.append(str(record.seq))

    protein_ed = editdistance.eval(fasta_sequence[0], fasta_sequence[1])

    print(protein_ed)

if __name__ == '__main__':
    main()
