#!/usr/bin/env python3
"""
@author: Timothy Baker
"""


freq_matrix = {'A':[0.01, 0.10, 0.97, 0.95, 0.50], 'C':[0.03, 0.05, 0.01, 0.01, 0.10], 'G':[0.95, 0.05, 0.01, 0.03, 0.10], 'T':[0.01, 0.80, 0.01, 0.01, 0.30]}


# I = 2 + sum(freq of the base and the position) * log2(frequency of the position and base)

import math


def calculate_information_content(frequency_matrix):

    position_to_calculate = int\
    (input("Enter position you want to calculate the Information Content --> "))

    position_to_calculate_index = position_to_calculate - 1
    sum_freq = 0

    for base in 'ACGT':
        prob = frequency_matrix[base][position_to_calculate_index]
        freq = prob*math.log2(prob)
        sum_freq += freq

    info_content = round(2 + sum_freq, 2)

    return info_content



def find_consensus_v3(frequency_matrix):
    """generates a consensus sequence given a frequency dictionary of lists"""

    consensus = ''
    dna_length = len(frequency_matrix['A'])

    for i in range(dna_length):  # loop over positions in string
        max_freq = -1            # holds the max freq. for this i
        max_freq_base = None     # holds the corresponding base

        for base in 'ACGT':
            if frequency_matrix[base][i] >= max_freq:
                max_freq = frequency_matrix[base][i]
                max_freq_base = base
            elif frequency_matrix[base][i] == max_freq:
                max_freq_base = '-' # more than one base as max

        consensus += max_freq_base  # add new base with max freq
    return consensus




def main():
    print(calculate_information_content(freq_matrix))
    print(find_consensus_v3(freq_matrix))



if __name__ == '__main__':
    main()
