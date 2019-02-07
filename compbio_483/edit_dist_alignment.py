#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 01-29-2019

edit_dist_alignment.py

Implementation Steps:
1) Read in FASTA file and parse sequences
2) Obtain global alignments of these sequences
3) Iterate through all alignments and calculate the minimum edit distance
    for each alignment using hamming distance
4) Store the alignment that has most minimum distance as the optimal alignment
    Could be multiple alignments, but unsure as to how to categorize them
    if all edit distances are the same.
Note:
    Weight the gap penalty more than mismatches.
"""

import sys
from Bio import SeqIO
from Bio import pairwise2



def hamming_distance(sequence_one, sequence_two):
    """ Calculates hamming distance between 2 strings
        Args:
            sequence_one (str) : sequence of ACGT
            sequence_two (str) : sequence of ACGT
        Returns:
            Hamming distance (int) : number of mismatches between 2 strings
        Raises:
            AssertionError : if sequence lengths are un equal
    """

    assert len(sequence_one) == len(sequence_two), "Sequence lengths are unequal"

    count = 0 #initialize the count to 0
    for n_1, n_2 in zip(sequence_one, sequence_two): #iterates through 2 sequences
        if n_1 != n_2: #checks equality against each nucleotide
            count += 1 # if not equal adds 1 to the count
    return count



def main():
    """ runs main script """

    # requires the fasta/txt file is in the same directory
    # must be in FASTA format
    fasta_text = sys.argv[1]

    # Biopython SeqIO Object
    fasta_list = list(SeqIO.parse(fasta_text, "fasta"))

    # Parse sequences to string objects
    fasta_sequence = [str(record.seq) for record in fasta_list]

    # generates global alignment with no match/mismatch, no gap penalty
    # pairwise2 does not produce a sorted list of tuples per source code
    alignments = pairwise2.align.globalms(fasta_sequence[0], fasta_sequence[1], 2, -1, -2, -2)

    # initializing the hamming distance using the first alignment
    # rosalind is using hamming distance as edit distance of equal length strings
    # alignment output tuples = (Algn1, Algn2, float, int, int)
    hamm_distance = hamming_distance(alignments[0][0], alignments[0][1])

    # initializing the optimal alignment with first aligned tuple
    optimal_alignment = alignments[0]

    # begin iterating through each alignment tuple generated from pairwise2
    # typically begins with start position
    for aligned_tuple in alignments:

        # calculate the hamming distance as a function of edit distance
        # strings are now equal length after alignment
        min_ham_edit_distance = hamming_distance(aligned_tuple[0], aligned_tuple[1])

        # if the distance is less than the initialized value,
        # reset the min_edit_distance and min alignment tuple
        # will not be replaced if an alignment has the same score
        if min_ham_edit_distance < hamm_distance:
            hamm_distance = min_ham_edit_distance
            optimal_alignment = aligned_tuple

    with open('output-edit.txt', 'w') as output:
        output.write(str(hamm_distance) + '\n' + \
        str(optimal_alignment[0]) + '\n' + \
        str(optimal_alignment[1]))


if __name__ == '__main__':
    main()
