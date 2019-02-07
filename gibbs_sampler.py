#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-04-2019

1E.py
GibbsSampler without Cromwell

Standard Input
    kmer number \n
    n iterations \n
    N times to run Gibbs function \n
    motif sequence 1 \n
    motif sequence n \n
Output > 1E-output.txt
    Best Gibbs Motifs: motif1...motifn \n
    Consensus Sequence: consensus motif \n
"""

import sys
import time
import random
from bisect import bisect

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



def construct_random_motif_list(dna_sequence_list, kmer_length):
    """ constructs a motif matrix of each first kmer from
        the dna_sequence_list
        Args:
            dna_sequence_list (lst) : list of nucleotide strings
            kmer_length (int) : kmer length
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
        i = random.randint(0, len(kmer_matrix[0]) - 1)
        best_kmer_motif.append(kmer_matrix[kmer_array[0]][i])

    return best_kmer_motif


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



def kmer_random_prob(prob_matrix, sequence, kmer_length):
    """ takes a probability dict, single string, and kmer_length, returns a kmer
        from the given sequence given an array of random probability distributions
        Args:
            prob_matrix (dict) : nucleotides are the sequence, values are an array
                of probabilities between 0 and 1.
            sequence (str) : sequence of nucleotides
            kmer_length (int) : desired kmer length
        Returns:
            kmer (str) : kmer string from kmer_seq at a randomly generated index

    """
    kmer_seqs = []
    total_probs = []
    init_tot = 0


    for i in range(len(sequence) - kmer_length + 1):

        kmer = sequence[i:i+kmer_length]

        # initializes the probability to float 1.0
        prob_tot = 1.0

        # multiples prob_tot against each random index, and updates the prob_tot
        # with the resulting product
        # if the probability of a specific nucleotide at some index is 0,
        # it zeroes out the prob_tot
        for ind, nuc in enumerate(kmer):
            prob_tot *= prob_matrix[nuc][ind]

        # sums the probability, does not have to equal 1.0
            init_tot += prob_tot

        kmer_seqs.append(kmer)
        total_probs.append(init_tot)

    # total_probs retains the sum of probabilities from each kmer
    # the following normalizes those sums back to between 0 and 1
    # from the those probabilities, a probability is randomly chosen
    try:
        norm_multiplier = 1.0 / sum(total_probs)
    except ZeroDivisionError:
        # if ZeroDivisionError, this will cause bisect to choose last element
        # from normal_total_probs, need to fix
        norm_multiplier = 0.0

    norm_total_probs = [i * norm_multiplier for i in total_probs]

    random_number = random.uniform(0, max(norm_total_probs))


    # utilizing binary search on norm_total_probs array and returns an index right of
    # whose value matches the random number chosen from an array of
    # normalized probabilities, then feeds the random index of kmer_seqs and
    # returns that kmer

    bin_index = bisect(norm_total_probs, random_number)
    return kmer_seqs[bin_index - 1]



def gibbs_sampler(dna_seq_list, kmer_length, n_times):
    """ takes an array of nucleotide sequences, desired kmer_length, and
        the number of times a random motif to be removed
        Args:
            dna_seq_list (lst) : list of strings
            kmer_length (int) : desired kmer length
            n_times (int) : number of times a motif is removed from the initial
                array
        Returns:
            best_motif_array (lst) : list of sequences of kmer length and that
                equal the length of the dna_seq_list

    """

    # the initial array of randomly selected kmers, 1 from each string
    # in dna_seq_list
    random_motif = construct_random_motif_list(dna_seq_list, kmer_length)

    # initializing the best_motif values with the initial selected random_motif
    # generates initial probability matrix, and scores for later use
    best_motif_array = random_motif
    best_motif_prob_matrix = build_profile(best_motif_array)
    best_motif_score = score_motifs(best_motif_prob_matrix, best_motif_array)


    # performs the sampling n_times
    # NOT to confuse with start points
    for j in range(n_times):

        # randomly selecting an integer up to the last index in dna_seq_list
        index_pos = random.randint(0, len(dna_seq_list)-1)

        # creates an array of the chosen motifs, while leaving the index_pos
        # selected motif out, does not change the size of random_motif during
        # iteration
        chosen_motifs = [motif for motif_index_pos, \
            motif in enumerate(random_motif) \
            if motif_index_pos != index_pos]

        # generates the probabilty matrix of the chosen motifs to use for the
        # kmer_random_prob function
        motif_prob_matrix = build_profile(chosen_motifs)

        # sets the new kmer randomly generated string to the
        # index_pos (ith) element of the random_motif
        random_motif[index_pos] = kmer_random_prob(motif_prob_matrix, \
                                                    dna_seq_list[index_pos], \
                                                    kmer_length)

        # print(random_motif) <- to check if the correct kmer at random
        # index position is being returned
        # computes the score for the random_motif array that contains the new
        # kmer_random_prob generated kmer string
        new_score = score_motifs(motif_prob_matrix, random_motif)

        # tests whether the score is lower than the initial and subsequent
        # best_motif_scores
        if new_score < best_motif_score:
            best_motif_array = random_motif
            best_motif_score = new_score

    return (best_motif_score, best_motif_array)





def main():
    """ takes standard input, runs main script, no error handling for bad input """

    lines = sys.stdin.read().splitlines()

    kmer_length = int(lines[0])
    n_times = int(lines[1])
    n_times_to_run_gibbs = int(lines[2])

    dna_sequence_list = lines[3:]

    best_gibbs_motifs = gibbs_sampler(dna_sequence_list, kmer_length, n_times)

    for i in range(n_times_to_run_gibbs-1):
        new_best_motifs = gibbs_sampler(dna_sequence_list, kmer_length, n_times)
        if new_best_motifs[0] < best_gibbs_motifs[0]:
            best_gibbs_motifs = new_best_motifs

    consensus_seq = find_consensus_sequence(build_profile(best_gibbs_motifs[1]))

    with open('1E-output.txt', 'w') as output:
        output.write("Best Gibbs Motifs: " + \
        " ".join(best_gibbs_motifs[1]) + \
        '\n' + "Consensus Sequence: " + \
        str(consensus_seq) + '\n')



if __name__ == '__main__':
    START_TIME = time.time()
    main()
    print(("--- %s seconds ---" % (time.time() - START_TIME)))
