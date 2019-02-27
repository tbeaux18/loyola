#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-25-2019

4B.py

Args to run file:
    -g/--gene_file : the supplement table csv file
    -d/--dog_file : FASTA format text/fasta file

Required:
    Supplemental_Table.csv that contains genotype ids with their sequences
    in [B/B] format at the specific loci

Input:
    FASTA fle containing dog contigs
Output:
    text file with ID-Allele that the sequence matches.

tested neg and positive cases, and partial sequence match cases
"""

import csv
import re
import argparse
from Bio import SeqIO

def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='Searches clinical file for concordant variants.\n '
    )
    parser.add_argument('-g', '--gene_file', help='Path to the gene Supplemental_Table.csv')
    parser.add_argument('-d', '--dog_file', help='Path to dog sequence file')

    return parser.parse_args()


class BurrowsWheelerTransform:
    """ Burrows Wheeler Transform Object for pattern matching. Takes a sequence
        and compresses to a burrows wheeler transformation, and a pattern to detect.

        Args:
            sequence (str) : DNA string object
            pattern (str) : DNA string object

        Attributes:
            sequence (str) : DNA string object from args
            pattern (str) : DNA string object from args
            bwt_seq (str) : initialized empty string object
            last_to_first (array) : initialized indexed array
        Methods:
            bwt() (instance) : rotates the sequence and returns the BWT sequence
            return_bwt() (instance) : returns the BWT sequence
            map_seqs() (instance) : last to first mapping to form index array
            return_map_idx_array() (instance) : returns the last to first map array
            bw_matching() (instance) : does the pattern matching against the
                                    compressed BWT sequence attribute


     """

    def __init__(self, sequence, pattern):
        self.sequence = sequence
        self.pattern = pattern
        self.bwt_seq = ""
        self.last_to_first = [None]


    def bwt(self):
        """ forms the BWT of the input sequence and returns it to map_seqs method """

        rotated_list = []

        for i in range(len(self.sequence)):
            rotated_list.append(self.sequence[i:]+self.sequence[:i])

        self.bwt_seq = ''.join([rot[-1] for rot in sorted(rotated_list)])

        return self.bwt_seq

    def return_bwt(self):
        """ returns the BWT sequence from the bwt() instance method """

        return self.bwt_seq()


    def map_seqs(self):
        """ performs the FM index mapping on the BWT sequence and returns
            the instantiated last_to_first array
        """
        self.bwt_seq = self.bwt()

        bwt_length = len(self.bwt_seq)

        # tuple each character with its index for ranking
        last_rank = [(self.bwt_seq[i], i) for i in range(bwt_length)]

        # sorting the last_rank yields the first rank, breaking ties
        # between same characters
        first_rank = sorted(last_rank)

        # creats an array of the ranked index
        first_to_last = [idx for (sym, idx) in first_rank]

        # initializing the last_to_first array with None objects for indexing
        # by multipling against the length of the BWT sequence.
        self.last_to_first = [None] * bwt_length

        # places rank's index value from first to last array into the
        # last to first array's index position that matches the rank value
        for firast in enumerate(first_to_last):
            self.last_to_first[firast[1]] = firast[0]

        return self.last_to_first

    def return_map_idx_array(self):
        """ returns the FM indexing array """

        return self.map_seqs()

    def bw_matching(self):
        """ performs the pattern matching against the BWT sequence using
            the last to first array for reconstruction
        """

        self.last_to_first = self.map_seqs()

        # initializing top pointer to 0
        top_pointer = 0

        # indexing the bottom pointer to the last index for the BWT length
        bottom_pointer = len(self.bwt_seq) - 1


        while top_pointer <= bottom_pointer:
            # will trigger false once the last character is ran
            if self.pattern:

                # grabs the last character of the pattern
                last_char = self.pattern[-1]
                # resets the pattern without the last character
                self.pattern = self.pattern[:-1]

                last2first_range = self.bwt_seq[top_pointer:bottom_pointer+1]

                # checking if the last character is in the BWT
                if last_char in last2first_range:

                    # locating the top pointer index if the char is present
                    top_index = self.bwt_seq.find(last_char, top_pointer, bottom_pointer+1)
                    # locating the bottom pointer index if the char is present
                    bottom_index = self.bwt_seq.rfind(last_char, top_pointer, bottom_pointer+1)

                    # resetting the top_pointer to the new top index
                    top_pointer = self.last_to_first[top_index]
                    # resetting the bottom pointer to the new bottom index
                    bottom_pointer = self.last_to_first[bottom_index]
                else:
                    return 0
            else:
                return bottom_pointer - top_pointer + 1


def handle_alleles(gene_dict):
    """ handles the Supplemental_Table.csv that holds the SNP IDs with their
        sequences. Sequences have one instance of [A/B] designating the
        nucleotide at that specific loci.
        Args:
            gene_dict (dct) : initialized from reading in CSV
        Returns:
            new_gene_dict (dct) : new dict containing the SNPID-A and SNPID-B
            alleles with their respective sequences for pattern matching
    """

    # regexing to split each sequence at the [A/B] location
    allele_regex = re.compile(r'\[(.*?)\]')

    new_gene_dict = {}

    for key, value in gene_dict.items():
        regex_list = allele_regex.split(value)

        # first part of the regex'd split sequence
        first = regex_list[0]
        # the allele site
        allele = regex_list[1].split('/')
        # the last part of the regex'd split sequence
        last = regex_list[2]

        # concatenate to form new string
        allele_a = first + allele[0] + last
        allele_b = first + allele[1] + last

        # add to new dict with each SNPID-A and SNPID-B and their sequences
        new_gene_dict[key+'-'+allele[0]] = allele_a
        new_gene_dict[key+'-'+allele[1]] = allele_b

    return new_gene_dict


def main():
    """ runs main script """

    args = arg_parser()

    gene_file = args.gene_file
    dog_file = args.dog_file

    dog_sequence_fasta = []

    for record in SeqIO.parse(dog_file, "fasta"):
        dog_sequence_fasta.append(str(record.seq))

    gene_dict = {}

    # opens the gene_file argument and performs handling
    with open(gene_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            gene_dict[row['SNP name']] = row['sequence']

    parsed_gene_dict = handle_alleles(gene_dict)

    genotype_hits = set()

    # iterates through the parsed gene dictionary and instantiates a
    # BWT object, the BWT'd sequence is the dog sequence and the patterns
    # to match are each genotype sequence

    for dog_sequence in dog_sequence_fasta:
        for key, values in parsed_gene_dict.items():

            bwm_obj = BurrowsWheelerTransform(dog_sequence, values)

            if bwm_obj.bw_matching():
                genotype_hits.add(key)

    with open('genotype-output.txt', 'w') as output:
        for geno_id in genotype_hits:
            output.write(geno_id+'\n')

if __name__ == '__main__':
    main()
