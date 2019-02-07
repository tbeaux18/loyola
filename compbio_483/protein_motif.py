t#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 10/19/2018

Finding a Protein Motif Rosalind Problem 14

"""

# missing regex pattern match

import sys
import re
from Bio import ExPASy
from Bio import SwissProt


def grab_motif_location(protein_id):
    """ Prints the locations of N-glycosylation motif from each record's
        protein sequence if present.
        TODO: Fix return formatting instead of printing to console
        Args:
            protein_id (str) : UniProt Accession ID
        Returns:
            None
    """

    output = [] #initializes record list

    # Pre-compiling regex pattern to match
    # Lookahead regex pattern to catch overlapping matches
    glyc_pattern = re.compile(r'(?=(N[^P][ST][^P]))')

    try:
        # passing the accession_id to grab the sequence record
        with ExPASy.get_sprot_raw(protein_id) as handle: #reads in the protein id
            record = SwissProt.read(handle)

            # locate pattern match within the record sequence
            regex = glyc_pattern.finditer(record.sequence)

            # iterates through each match group and appends start index position
        for match in regex:
            output.append(match.start(1) + 1)


        if len(output) >= 1:
            #return "{}\n{}".format(protein_id, " ".join('%d'*(len(output)) % tuple(output)))
            print(protein_id)
            print(" ".join(str(i) for i in output))


    except ValueError:
        print("WARNING: Accession %s not found" % protein_id)





def main():
    """ obtains list of ids from std input and prints each id with motif locations """

    uniprot_ids = [x.rstrip() for x in sys.stdin.readlines()]

    for accession_id in uniprot_ids:
        grab_motif_location(accession_id)



if __name__ == '__main__':
    main()
