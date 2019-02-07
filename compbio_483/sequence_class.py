#!/usr/bin/env python3

<<<<<<< HEAD
TEST_SEQ = 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'


class Sequence():

    def __init__(self, sequence):
        self.sequence = sequence


    def count_nucleotide(self):
        """ takes a sequence object and counts ACGT """

        a_count, c_count, g_count, t_count = 0, 0, 0, 0

        for nucleotide in self.sequence.upper():
=======

# Rosalind challenge #1
# @author: Timothy Baker
# @version: 22/08/2018

# Given: A DNA string s of length at most 1000 nt.
# Return: Four integers (separated by spaces) counting the
#        respective number of times that the symbols 'A', 'C', 'G', and 'T'
#        occur in s.

import sys


class Sequence:
    """ DNA Sequence object

        Attributes:
        sequence: a dna sequence string

    """

    def __init__(self, nt_sequence):
        """ inits a string sequence """

        self.nt_sequence = nt_sequence.upper()


    def dna_to_rna(self):
        """ converts dna to rna """

        rna = self.nt_sequence.replace('T', 'U')
        return rna


    def generate_complement(self):
        """ generates the complementary strand """

        for nucleotide in self.nt_sequence:
            if nucleotide == 'A':
                self.nt_sequence.replace('A', 'T')
            elif nucleotide == 'C':
                self.nt_sequence.replace('C', 'G')
            elif nucleotide == 'G':
                self.nt_sequence.replace('G', 'C')
            else:
                self.nt_sequence.replace('T', 'A')

        return self.nt_sequence



    def count_nt_sequence(self):
        """ Performs nt counting """

        a_count, c_count, g_count, t_count = 0, 0, 0, 0

        for nucleotide in self.nt_sequence:
>>>>>>> 89b4129ffcc8d394f850283f63ae3dc0749372e6
            if nucleotide == 'A':
                a_count += 1
            elif nucleotide == 'C':
                c_count += 1
            elif nucleotide == 'G':
                g_count += 1
            else:
                t_count += 1
<<<<<<< HEAD

        return (a_count, c_count, g_count, t_count)


    def transcribe_dna_to_rna(self):

        rna_seq = [n for n in self.sequence.upper() if n ==''
        for nuc in self.sequence.upper():
            if nuc == ''

def main():
    test_seq_object = Sequence(TEST_SEQ)



=======
        return a_count, c_count, g_count, t_count



def main():
    """ Parses standard input file and creates an instance of a
        type DnaSequence object. Performs a nt count on the string.
        Outputs to a standard output file.
    """
    with open('rosalind_dna.txt', 'r') as raw_seq:
        for line in raw_seq:
            raw_sequence = line.rstrip()
    dnaseq_object = Sequence(raw_sequence)
    #count_dnaseq_object = dnaseq_object.count_nt_sequence()
    #rna_object = dnaseq_object.dna_to_rna()
    comp_test = dnaseq_object.generate_complement()
    print(comp_test)
>>>>>>> 89b4129ffcc8d394f850283f63ae3dc0749372e6

if __name__ == '__main__':
    main()
