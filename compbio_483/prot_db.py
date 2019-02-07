#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 10/16/2018

Introduction to Protein Databases Rosalind Problem 12
"""

import sys
from Bio import ExPASy
from Bio import SwissProt


def grab_prot_bioprocess(protein_id):
    """ Grabs the biological processes from uniprotkb for a given record
        Args:
            protein_id (str) : std input into a string; will return
            HTTP 404 Error if id not found in database
        Returns:
            biological processes from the gene ontology section of unitkb record
    """

    with ExPASy.get_sprot_raw(protein_id) as handle: #reads in the protein id
        record = SwissProt.read(handle) #stores handle as a string

    gene_ont_list = [] #initialize a list

    # iterates through tuple pairs in cross_references
    for tuple_pair in record.cross_references:
        # checks if zero index is GO for Gene Ontology
        # checks if second index is P for Biological Process
        if tuple_pair[0] == 'GO' and tuple_pair[2].startswith('P'):
            # append the biological process to gene ont list
            gene_ont_list.append(tuple_pair[2][2:])

    return "\n".join(x for x in gene_ont_list) #returns in newline format

def main():
    """ Takes protein id from std input and prints the bio process to console """

    prot_id = sys.stdin.read().rstrip()
    print(grab_prot_bioprocess(prot_id))




if __name__ == '__main__':
    main()
