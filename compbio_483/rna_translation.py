#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 11/10/2018

Translating RNA into Protein Rosalind Problem 9

"""

import sys


def construct_aa_dict(file_name):
    """ Takes a file in the format CODON\tAMINO ACID ABBREV
    and converts this curated file to a dictionary. Not reproducible
    without file in this format.
        Args:
            file_name (txt) :
        Returns:
            dict : key = CODON; values = AMINO ACID ABBREV

    """
    with open(file_name, 'r') as aa_table:

        # takes a string from the file and converts it to a list of lists.
        aa_table_read = [",".join(x.split()) for x in aa_table.readlines()]

        aa_dict = {}

        for amino in aa_table_read:
            aa_dict[amino.split(',')[0]] = amino.split(',')[1]

        return aa_dict


def protein_translation(rna_sequence, amno_dict, kmer_length=3):
    """
        Args:
            rna_sequence (str) : Must contain Uracil
            kmer_length (int) : set to 3 as default
            amno_dict (dict) : set to amino dict from construct_aa_dict function
        Returns:
            translated protein sequence
    """

    codon_list = [] # initialize the codon list

    protein_sequence = [] # initialized the protein sequence list


    for i in range(0, len(rna_sequence), 3): #iterates through the str sequence
        codon_list.append(rna_sequence[i:i+kmer_length]) #generates the codon and appends to a list


    for codon in codon_list: #iterates over each codon
    #adds amino acids according to the amino acid table dict that does not equal Stop
    # Does not work if a Stop codon is in middle of the sequence
        if amno_dict[codon] != 'Stop':
            protein_sequence.append(amno_dict[codon])

    return "".join(protein_sequence) #joins the protein sequence list into a string



def main():
    """ Constructs the amino acid table, takes in std input format of rna rna_sequence
    and prints the rna to protein translation to the console
    """

    aa_table_name = 'dna_codon_table.txt'
    amino_dict = construct_aa_dict(aa_table_name)

    #rna_line = sys.stdin.read().strip()
    dna = 'CAGGCACCAGGGACATTATTTACGTTTCCGTCACTCACG'

    print(protein_translation(dna, amno_dict=amino_dict))



if __name__ == '__main__':
    main()
