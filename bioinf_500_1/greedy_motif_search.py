#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 01-30-2019

1D.py

Input:
    kmer number \n
    motif sequence 1 \n
    motif sequence n \n
Output:
    Best Motifs: motif1...motifn \n
    Consensus Sequence: consensus motif \n
"""

import sys
import time
import functools
import operator

def generate_kmer(dna_sequence, kmer_length):
    """ Takes a DNA sequence and generates all kmer substrings of the sequence
        Args:
            dna_sequence (str) : must be ACGT; no RNA transcripts
            kmer_length (int) : default is set to 3 for codon, but can be changed
        Returns:
            a list of all possible kmers
    """

    kmer_list = [] # initializes the kmer list; order is important

    # iterates through the sequence creating a kmer length sliding window
    if isinstance(dna_sequence, list):
        for seq in dna_sequence:
            for i in range(len(seq) - kmer_length + 1):
                # appends each kmer to the kmer list
                kmer_list.append(seq[i:i+kmer_length])


    if isinstance(dna_sequence, str):
        for i in range(len(dna_sequence) - kmer_length + 1):
            # appends each kmer to the kmer list
            kmer_list.append(dna_sequence[i:i+kmer_length])

    return kmer_list



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



def construct_best_motif_list(dna_sequence_list, kmer_length):
    """ constructs a motif matrix of each first kmer from
        the dna_sequence_list
        Args:
            dna_sequence_list (lst) : list of nucleotide strings
            kmer_length (int) : kmer length
            best (bool) : selects the first kmer from each kmer array
        Returns:
            kmer_motif_matrix (txk matrix) : first kmer from each string
                forming a t x k matrix
    """

    kmer_matrix = []
    best_kmer_motif = []

    for seq in dna_sequence_list:
        kmer_matrix.append(generate_kmer(seq, kmer_length))

    # this selects the first kmer from each array
    for kmer_array in enumerate(kmer_matrix):
        best_kmer_motif.append(kmer_matrix[kmer_array[0]][0])

    return best_kmer_motif


def build_profile(motif_list):
    """ generates a dictionary with frequency values at each position
        Args:
            motif_list (lst) : list of sequences of equal length
        Returns:
            profile (dict) : returns a dict with probability values
                keys are nucleotides, should only ever be 4 keys

    """

    k = len(motif_list[0])

    profile = {'A':[0]*k, 'C':[0]*k, 'G':[0]*k, 'T':[0]*k}

    div = float(len(motif_list))

    for i in range(k):
        for motif in motif_list:
            profile[motif[i]][i] += 1
        for key in profile:
            profile[key][i] /= div

    return profile



def prod(factors):
    """ reduces the multipled probability of the probability matrix to 1 or 0"""
    return functools.reduce(operator.mul, factors, 1)



def return_kmer_most_prob(kmer_list, prob_matrix):
    """ takes a list of kmers, and a probability matrix, and returns
        a dict of each kmer, and uses the prod function to label the
        kmer with the highest probability (1)
        Args:
            kmer_list (lst) : list of kmers of equal length
            prob_matrix (dict) : dict from build_profile function with the
            4 nucleotides as they keys and a list of probabilities to be
            multiplied
        Returns:
            dict with kmers as keys, and the value 1 or 0 depending on the
            highest multipled probability.
    """
    return {
        kmer: prod(
            prob_matrix[base][i]
            for i, base in enumerate(kmer)
        )
        for kmer in kmer_list
    }



def build_profile_probable(dna_sequence, kmer_length, prob_matrix):
    """ takes in a sequence, generates the kmers of that sequence, and
        a probabilitiy matrix.
        Utilizes return_kmer_most_prob that multiples the probabilities
        and reduces to a single value either 0 or 1.
        Args:
            dna_sequence (str) : sequence
            kmer_length (int) : kmer_length
            prob_matrix (dict of lsts) : typically from build_profile function
            4 nucleotides with a list of probabilities
        Returns:
            most_prob_motif (str) : the key (kmer) with 1 as the value

    """

    kmer_list = generate_kmer(dna_sequence, kmer_length)

    most_prob_kmer = return_kmer_most_prob(kmer_list, prob_matrix)

    most_prob_motif = max(most_prob_kmer.items(), key=operator.itemgetter(1))[0]

    return most_prob_motif



def find_consensus_sequence(probability_matrix):
    """generates a consensus sequence given a frequency dictionary of lists
        Args:
            frequency_matrix (dict of lists) : probability matrix

        Returns:
            consensus (str) : baesd on the probability matrix, constructs a
            consensus string
    """

    consensus = ''
    dna_length = len(probability_matrix['A'])

    for i in range(dna_length):  # loop over positions in string
        max_freq = -1            # holds the max freq. for this i
        max_freq_base = None     # holds the corresponding base

        for base in 'ACGT':
            if probability_matrix[base][i] >= max_freq:
                max_freq = probability_matrix[base][i]
                max_freq_base = base
            elif probability_matrix[base][i] == max_freq:
                max_freq_base = '-' # more than one base as max

        consensus += max_freq_base  # add new base with max freq
    return consensus


def score_motifs(prob_matrix, motif_list):
    """ takes in probability matrix, and a motif list, and returns a score
        calculates the hamming distance between each motif and the consensus
        sequence.
        Args:
            prob_matrix (dict of lists) : prob matrix with 4 nucleotide keys
            and a list of probabilities
            motif_list (lst) : list of motifs of equal length
        Returns:
            score (int) : the summation of the hamming_distance between the consensus
            sequence and each motif

    """
    score = 0

    consensus_sequence = find_consensus_sequence(prob_matrix)

    for motif in motif_list:
        score += hamming_distance(consensus_sequence, motif)

    return score



def greedy_motif_search(dna_seq_list, kmer_length):
    """ takes in a list of dna_sequences, and a desired kmer_length, returns
        an array of kmers with the lowest hamming distance score between
        a constructed consensus sequence and motifs
        Args:
            dna_seq_list (lst) : list of sequences
            kmer_length (int) : desired kmer length, cannot be longer
            than the sequences entered
        Returns:
            best_motif_array (lst) : an array of the lowest scoring motifs

    """

    t_length = len(dna_seq_list)

    # initializing best motifs by taking the first kmer from each sequence
    # and forming a probability matrix, and score
    best_motif_array = construct_best_motif_list(dna_seq_list, kmer_length)
    best_motif_prob_matrix = build_profile(best_motif_array)
    best_motif_score = score_motifs(best_motif_prob_matrix, best_motif_array)

    # iterating through each kmer of the first sequence in the list
    for kmer_motif in generate_kmer(dna_seq_list[0], kmer_length):

        # initializing each kmer from the first sequence to an array
        kmer_motif_array = [kmer_motif]

        # iterating through the second to the last sequence
        for i in range(1, t_length):

            # first pass generates a probability matrix using a single kmer from
            # the first element in list from the generate_kmer function
            # at each iterating pass, the kmer most probable is added to the
            # kmer_motif_array, updating the probability matrix to be more accurate
            prob_matrix = build_profile(kmer_motif_array)

            most_prob_kmer_motif = build_profile_probable(dna_seq_list[i], kmer_length, prob_matrix)

            kmer_motif_array.append(most_prob_kmer_motif)

        # constructs a score for each kmer_motif_array after reaching the t_length
        # and compares the score against best_motif_score, if less than, the best
        # scoring motif array is updated with the lowest score.
        new_score = score_motifs(prob_matrix, kmer_motif_array)
        if new_score < best_motif_score:
            best_motif_array = kmer_motif_array
            best_motif_score = new_score


    return best_motif_array



def main():
    """ runs mains script """

    lines = sys.stdin.read().splitlines()

    kmer_length = int(lines[0])

    motif_sequence_list = lines[1:]

    motif_array = greedy_motif_search(motif_sequence_list, kmer_length)

    consensus_seq = find_consensus_sequence(build_profile(motif_array))

    with open('1D-output.txt', 'w') as output:
        output.write("Best motifs: " + \
        " ".join(motif_array) + \
        '\n' + \
        "Consensus Sequence: " + \
        str(consensus_seq) + '\n')

if __name__ == '__main__':
    START_TIME = time.time()
    main()
    print(("--- %s seconds ---" % (time.time() - START_TIME)))
