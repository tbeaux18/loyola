#!/usr/bin/env python3
"""
@author: Timothy Baker

seq_gen.py

Random functions that can be called to generate random DNA sequences of x length

"""


def build_dna_string(length):
    """ generates a single x length nucleotide sequence on 'ACGT' at random
        Args:
            length (int) : length of sequence
        Returns:
            nucleotide_sequence (str) : sequence of 'ACGT'
    """
    from random import choice

    dna_sequence = ""

    for count in range(length):
        dna_sequence += choice("ACGT")

    return dna_sequence






def replace_string(sequence, motif, index, nofail=False):
    """ takes input string, motif sequence, and index and creates a new
        string
        Args:
            sequence (str) : any string
            motif (str) : any new string needed to be added
            index (int) : cannot be greater than the length of sequence
            nofail (bool) : default False
        Return:
            concatenated sequence (str)
    """
    # raise an error if index is outside of the string
    if not nofail and index not in range(len(sequence)):
        raise ValueError("index outside given string")

    # if not erroring, but the index is still not in the correct range..
    if index < 0:  # add it to the beginning
        return motif + sequence
    if index > len(sequence):  # add it to the end
        return sequence + motif

    # insert the new string between "slices" of the original
    return sequence[:index] + motif + sequence[index + 1:]





def implant_motif(dna_sequence_list, motif_length):
    """ takes a list of sequences, and a motif length and inputs the SAME
        motif into random index positions
        Args:
            dna_sequence_list (lst) : list of strings
            motif_length (int) : any integer length
        Returns:
            motif_seq_list (lst) : list of same sequences with motifs inserted
            into random positions
    """
    from random import randint

    motif_seq_list = []

    # only 1 motif per sequence list
    # motif prints to console to alert user of the sequence
    motif_sequence = build_dna_string(motif_length)
    print(motif_sequence)


    for seq in dna_sequence_list:
        seq_length = len(seq) - 1
        rand_index = randint(1, seq_length)
        motif_seq_list.append(replace_string(seq, motif_sequence, rand_index))

    return motif_seq_list



def main():
    """ Example code below to use functions"""

    # If just creating random DNA sequences of x length the below
    # list comprehension will generate a list of sequences.
    sequence_list = [build_dna_string(20) for i in range(1, 3)]
    # If wanting to create sequences with single motif implanted into
    # each sequence in the sequence_list.
    motif_sequence_list = implant_motif(sequence_list, 10)

    print(sequence_list)

if __name__ == '__main__':
    main()
