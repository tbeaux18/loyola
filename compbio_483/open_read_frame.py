#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 10/16/2018

Open Reading Frames Rosalind Problem 13

"""

import sys
import reverse_complement as rc
import rna_translation as rt




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



def translate_dna(kmer_list, start_position, amino_dict):
    """ takes a list of all kmers, and the start position of the start codon,
        and returns a joined string of the start codon and every third codon
        Args:
            kmer_list (list) :
            start_position (int) :
            amino_dict (dict) :
        Returns:
            protein_seq (list) : only every third codon
    """

    # initializes the protein_seq list
    protein_seq = []

    # iterates through the kmer_list
    for codon in enumerate(kmer_list[start_position:]):

        # checks if the index of the codon is divisible by 3 and only adds every
        # third codon since each possible kmer was generated
        if codon[0] % 3 == 0:
            protein_codon = amino_dict[codon[1]] # sets the protein from the amino acid dict
            if protein_codon == 'Stop':
                return protein_seq # returns the protein without the Stop protein
            else:
                # adds to protein until Stop is reached
                protein_seq.append(protein_codon)

    return protein_seq






def open_frame_translation(sequence, amino_dict):
    """ Translates a DNA sequence by translating all open frames to protein.
        Start Codons: ATG (M)
        Stop Codons: TAA, TAG, TGA (Stop)
        Args:
            sequence (str) : DNA sequence, must be ACGT
            amino_dict (dict) : amino acid dna codon table from diff. function
        Returns:
            translated proteins for all open frames (list)

    """

    # initializes the open_protein_seq list
    open_protein_seq = []

    # generates all 3-mers for the given sequence
    all_kmer_list = generate_kmer(sequence)

    start_position = all_kmer_list.index('ATG') # initialize start position

    try:
        # while statement, index method returns a ValueError when it gets to
        # the end of the sequence
        while start_position != ValueError:

            # creates the protein product from the translate function
            protein_product = translate_dna(all_kmer_list, start_position, amino_dict)

            # deletes the used start codon, shouldn't effect the sequence
            del all_kmer_list[start_position]

            # reinitializes a start position after the previous one is deleted
            start_position = all_kmer_list.index('ATG') # re-starts start position

            # Adds each protein product to the open_protein_seq list
            if protein_product is not ValueError:
                open_protein_seq.append(protein_product)


    except ValueError: # raises ValueError and returns the final open_protein_seq
        return open_protein_seq


def handle_output(output_list_one, output_list_two):
    """ takes two lists and adds them to a set to deduplicate them
        Args:
            output_list_one (list)
            output_list_two (list)
        Returns:
            set of deduplicated proteins (set)
    """

    final_set = set() # initialize the set

    for pro_one in output_list_one: # iterates through the list
        final_set.add("".join(pro_one)) # adds to the set


    for pro_two in output_list_two: # iterates through the list
        final_set.add("".join(pro_two)) # adds to the set

    return "\n".join(final_set) #returns newline format



def main():
    """ Creates amino acid dictionary, parses standard input fasta record,
    and prints the translated open frames to console
    """

    # generates the amino acid dictionary from specific dna codon table
    aa_dict = rt.construct_aa_dict('dna_codon_table.txt')

    # standard input fasta record, splits by > title and sequence
    fasta_record = [x.rstrip() for x in sys.stdin]

    # concatenated sequence
    fasta_sequence = "".join(x for x in fasta_record if not x.startswith('>'))

    # runs the translation function
    protein_coding_strand = open_frame_translation(fasta_sequence, aa_dict)
    protein_complement_strand = open_frame_translation(\
    rc.generate_reverse_complement(fasta_sequence), aa_dict)

    # prints the output to the console in string newline format
    with open('orf_output.txt', 'w') as orf_output:
        orf_output.write(handle_output(protein_coding_strand, protein_complement_strand))
    print(handle_output(protein_coding_strand, protein_complement_strand))

if __name__ == '__main__':
    main()
