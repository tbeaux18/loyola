#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 01-29-2019

1C.py

Standard input
    kmer number \n
    d mismatch number \n
    motif sequence 1 \n
    motif sequence n \n

Output format
    Motif sequences detected: motif1 ... motifn
"""

import sys
import time


def generate_kmer(dna_sequence, kmer_length, d_mismatch, neighbors=False):
    """ Takes a DNA sequence and generates all kmer substrings of the sequence
        Args:
            dna_sequence (str) : must be ACGT; no RNA transcripts
            kmer_length (int) : default is set to 3 for codon, but can be changed
        Returns:
            a list of all possible kmers
    """

    kmer_list = [] # initializes the kmer list; order is important

    list_of_neighbors_from_kmers = []

    # iterates through the sequence creating a kmer length sliding window
    if isinstance(dna_sequence, list):
        for seq in dna_sequence:
            for i in range(len(seq) - kmer_length + 1):
                # appends each kmer to the kmer list
                kmer_list.append(seq[i:i+kmer_length])
        return kmer_list

    if isinstance(dna_sequence, str):
        if neighbors is False:
            for i in range(len(dna_sequence) - kmer_length + 1):
                # appends each kmer to the kmer list
                kmer_list.append(dna_sequence[i:i+kmer_length])
            return kmer_list

        if neighbors is True:
            for i in range(len(dna_sequence) - kmer_length + 1): #
                list_of_neighbors_from_kmers.append(\
                build_neighborhood(dna_sequence[i:i+kmer_length], d_mismatch))

            return list_of_neighbors_from_kmers



def build_neighborhood(kmer_pattern, d_mismatch):
    """ recursively returns list of all kmer permutions for a given kmer_pattern
        with at most d-mismatches
        Args:
            kmer_pattern (str) : any length str
            d_mismatch (int) : any number of integers that do not exceed kmer length
        Returns:
            set of all kmer permutations
    """

    if d_mismatch == 0:
        return kmer_pattern

    elif len(kmer_pattern) == 1:
        return ['A', 'C', 'G', 'T']

    neighborhood = []

    suffix_neighbor = build_neighborhood(kmer_pattern[1:], d_mismatch)

    for kmer_suffix in suffix_neighbor:
        if hamming_distance(kmer_pattern[1:], kmer_suffix) < d_mismatch:
            for nucleotide in ['A', 'C', 'G', 'T']:
                neighborhood.append(''.join([nucleotide, kmer_suffix]))
        else:
            neighborhood.append(''.join([kmer_pattern[0], kmer_suffix]))

    return neighborhood



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


def iter_flatten(iterable):
    """ flattens a list of lists """

    iterd = iter(iterable)

    for ele in iterd:
        if isinstance(ele, (list, tuple)):
            for flt in iter_flatten(ele):
                yield flt
        else:
            yield ele


def motif_enumeration(sequence_list, kmer_length, d_mismatch):
    """ returns the intersection of all generated kmer neighborhood sets
        with at most d-mismatches produced by build neighborhood function
        Args:
            sequence_list (lst) : list of sequences with motif implanted
            kmer_length (int) : length of the kmer, typically should equal motif length
            d_mismatch (int) : number of allowed mismatches in order to produce
                                kmer permutations
            Returns:
                intersection of all sets else prints nothing found
    """

    k_mer_list = []

    new_kmer_list_w_neighbors = []
    for seq_tup in enumerate(sequence_list):
        k_mer_list.append(generate_kmer(sequence_list[seq_tup[0]], \
                                        kmer_length, \
                                        d_mismatch, \
                                        neighbors=True))


    for i in range(len(k_mer_list)):
        new_kmer_list_w_neighbors.append([lst for lst in iter_flatten(k_mer_list[i])])
    #print(new_kmer_list_w_neighbors)
    intersection = set.intersection(*map(set, new_kmer_list_w_neighbors))

    return intersection


 #  __  __          _____ _   _
 # |  \/  |   /\   |_   _| \ | |
 # | \  / |  /  \    | | |  \| |
 # | |\/| | / /\ \   | | | . ` |
 # | |  | |/ ____ \ _| |_| |\  |
 # |_|  |_/_/    \_|_____|_| \_|


def main():
    """ runs main script """


    lines = sys.stdin.read().splitlines()

    kmer_length = int(lines[0])
    dmismatch = int(lines[1])
    motif_sequence_list = lines[2:]

    motif_intersection = motif_enumeration(motif_sequence_list, kmer_length, dmismatch)

    pretty_new_intersection = ' '.join(str(item) for item in motif_intersection)

    with open('1C-output.txt', 'w') as output:
        output.write("Motif sequences detected: " + pretty_new_intersection + '\n')


if __name__ == '__main__':
    START_TIME = time.time()
    main()
    print(("--- %s seconds ---" % (time.time() - START_TIME)))
