#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-18-2019

assemble_shortest_contig.py

Dependencies:
    biopython

"""

import sys
from Bio import SeqIO

def assemble_overlap(sequence_list, contig=""):
    """ takes in a list of sequences and instantiates an empty string contig
        Args:
            sequence_list (lst) : list of strings
            contig (str) : empty
        Returns:
            recursively concatenates elements of the list together so long
            as sequence A's suffix matches sequence B's prefix
    """

    # if the list is 0, return the contig
    # occurs after last 2 elements are concatenated together
    if not sequence_list:
        return contig

    # initializes the contig to the first element in the list
    if not contig:
        contig = sequence_list.pop(0)
        return assemble_overlap(sequence_list, contig)

    else:
        # iterates through the elements
        for idx, seq in enumerate(sequence_list):
            length = len(seq)

            # begins to index prefix sequences half the length of the sequence
            for prefix in range(length // 2):
                diff = length - prefix

                # if the initial and subsequent contigs prefix starts with the sequence
                # suffix then it removes the ith element from the list and recursively
                # calls itself to concatenate the contig and iterating sequence together
                if contig.startswith(seq[prefix:]):
                    sequence_list.pop(idx)
                    return assemble_overlap(sequence_list, seq[:prefix] + contig)

                # else if the contig's suffix is the iterative sequence's prefix
                # removes the element and recursively calls itself to concatenate
                # the contig and iterating sequence together
                if contig.endswith(seq[:diff]):
                    sequence_list.pop(idx)
                    return assemble_overlap(sequence_list, contig + seq[diff:])

def main():
    """ runs main script """

    # requires the fasta/txt file is in the same directory
    fasta_text = sys.argv[1]

    fasta_list = list(SeqIO.parse(fasta_text, "fasta"))

    sequence_list = [str(record.seq) for record in fasta_list]

    result = assemble_overlap(sequence_list)

    with open('overlap-output.txt', 'w') as output:
        output.write(str(result))
if __name__ == '__main__':
    main()
