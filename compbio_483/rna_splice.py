#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 01-22-2019

rna_splice.py

Input:


"""
import sys
from Bio import SeqIO
import rna_translation as rnat


def split_string(string, intron):
    """ to split string on intron """

    # Split the string based on space delimiter
    list_string = string.split(intron)

    return list_string

def join_string(list_string):
    """ joins new exon post intronic split """

    # Join the string based on no delimiter
    string = ''.join(list_string)

    return string


def cut_introns(fasta_list):
    """ takes fasta_list and returns the concatenated exons sequenceself.
        Ugly fix, exon_list must remain the same to overwrite existing value
        Args:
            fasta_list (SeqIO)
        Returns:
            exon_list (str)

    """
    total_sequences = []

    for record in fasta_list:
        total_sequences.append(str(record.seq))

    exon_list = total_sequences[0]
    intron_list = total_sequences[1:]

    for intron in intron_list:
        exon_list = join_string(split_string(exon_list, intron))


    return exon_list



def main():
    """ runs main script """
    fasta_text = sys.argv[1]

    fasta_list = list(SeqIO.parse(fasta_text, "fasta"))

    exon_sequence = cut_introns(fasta_list)

    aa_table_name = 'dna_codon_table.txt'
    amino_dict = rnat.construct_aa_dict(aa_table_name)

    print(rnat.protein_translation(exon_sequence, amno_dict=amino_dict))

if __name__ == '__main__':
    main()
