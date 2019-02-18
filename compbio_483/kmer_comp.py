#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 01-22-2019

kmer_compy.py

"""

import sys
import collections
from Bio import SeqIO

def generate_kmer(dna_sequence, kmer_length=3):
    """ Takes a DNA sequence and generates all kmer substrings of the sequence
        Args:
            dna_sequence (str) : must be ACGT; no RNA transcripts
            kmer_length (int) : default is set to 3 for codon, but can be changed
        Returns:
            a list of all possible kmers
    """

    kmer_list = [] # initializes the kmer list; order is important

    # iterates through the sequence creating a kmer length sliding window
    for i in range(len(dna_sequence) - kmer_length + 1):
        # appends each kmer to the kmer list
        kmer_list.append(dna_sequence[i:i+kmer_length])

    return kmer_list


def main():
    """ runs main script """

    fasta_text = sys.argv[1]

    fasta_list = list(SeqIO.parse(fasta_text, "fasta"))

    fasta_sequence = []
    for record in fasta_list:
        fasta_sequence.append(str(record.seq))

    print("Generating kmers...")
    new_kmer_list = sorted(generate_kmer(fasta_sequence[0], kmer_length=4))

    print("Counting instances of each kmer")
    kmer_count_dict = collections.Counter(new_kmer_list)

    print("Sorting the keys...")
    kmer_order_dict = collections.OrderedDict(sorted(kmer_count_dict.items()))

    print("Writing to txt file...")
    with open('kmer-comp.txt', 'w') as output:
        for value in kmer_order_dict.values():
            output.write(str(value) + '\t')

if __name__ == '__main__':
    main()
