#!/usr/bin/env python3

"""
@author: Timothy Baker
@date: 01-24-2019

1A.py

Standard Input:
    Text file received as standard input
        Line 1) kmer length (k); any integer
        Line 2) number of occurrences (t); any integer
        Line 3) DNA sequence
Output:
    Text file showing Kmer: {Kmer Sequence} Count: {Kmer Count}
Exceptions:
    If there are no kmers to present meeting the number of occurence conditions
    line is printed to console stating such.
"""

import sys
import collections
import time


class Sequence:
    """ Sequence object that builds kmer list and counts and returns the number of
        kmers based on threshold input
    """

    def __init__(self, sequence, kmer_length, threshold):
        """ default constructor for sequence object
            Args:
                sequence (str): dna sequence
                kmer_length (int): desired kmer length (k)
                threshold (int): desired threshold for kmer count (t)
        """
        self.sequence = sequence
        self.kmer_length = kmer_length
        self.threshold = threshold


    def generate_kmer(self):
        """ generates list of kmers of length k """
        kmer_list = []

        for i in range(len(self.sequence) - self.kmer_length + 1):
            kmer_list.append(self.sequence[i:i+self.kmer_length])

        return kmer_list


    def count_kmer(self):
        """ creates a kmer count dict of all kmers """
        kmer_count_dict = collections.Counter(self.generate_kmer())

        return kmer_count_dict


    def filter_threshold(self):
        """ filters kmer count dict to a new dict based on the threshold """
        kmer_thres = {}

        for kmer_key, count_value in self.count_kmer().items():
            if count_value > self.threshold:
                kmer_thres[kmer_key] = count_value

        return kmer_thres


def main():
    """ takes std input, builds sequence object, writes kmer threshold dict to
        output txt file. results in same directory.
    """

    lines = sys.stdin.read().splitlines()

    if len(lines) > 3:
        kmer_length = lines[0]
        threshold_occurrence = lines[1]
        dna_sequence = ''.join(lines[2:])


    elif len(lines) == 3:
        kmer_length = lines[0]
        threshold_occurrence = lines[1]
        dna_sequence = lines[2]

    # builds sequence object
    seqobj = Sequence(dna_sequence, int(kmer_length), int(threshold_occurrence))


    # kmer count dict
    kmer_threshold = seqobj.filter_threshold()


    # handles unmet thresholds
    if kmer_threshold:
        # output to txt file
        print("Writing to output file...")
        with open('1A-output.txt', 'w') as output:
            for key, value in kmer_threshold.items():
                output.write('Kmer: {}\tCount: {}\n'.format(key, value))
    else:
        print("No {}-mers present in the sequence that meet the threshold value {}.".format(\
        kmer_length, threshold_occurrence))



if __name__ == '__main__':
    START_TIME = time.time()
    main()
    print(("--- %s seconds ---" % (time.time() - START_TIME)))
