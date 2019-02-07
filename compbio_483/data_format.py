#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 10/15/2018

Data Formats Rosalind Problem 11

"""

import sys
from Bio import Entrez
from Bio import SeqIO



def retrieve_minseq_fasta(genbank_id_list):
    """ retrieves fasta file with minimum sequence length out of list
        Args:
            genbank_id_list (int) : list of genbank ids
        Returns:
            FASTA formatted record that has minimum sequence length
    """
    Entrez.email = "tbeaux18@gmail.com" #credentialing

    #fetch the fasta files from the genbank_id_list
    handle = Entrez.efetch(db="nucleotide", id=genbank_id_list, rettype="fasta")

    #converts to a list of SeqRecord objects
    records_list = list(SeqIO.parse(handle, "fasta"))

    #creates lengths of each SeqRecord.seq object
    seq_length_list = [len(rec.seq) for rec in records_list]

    #finds the index with minimum length
    min_index = seq_length_list.index(min(seq_length_list))

    #returns fasta format of record with minimum seq length
    return records_list[min_index].format("fasta")




def main():
    """ Creates genbank list from standard input and prints the fasta file with
    minimum sequence length in fasta format
    """

    for line in sys.stdin:
        genbank_list = line.split()
    print(retrieve_minseq_fasta(genbank_list))


if __name__ == '__main__':
    main()
