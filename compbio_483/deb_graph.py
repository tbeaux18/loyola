#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-10-2019

deb_graph.py

Constructing a De Bruijn Graph

"""

import sys
from Bio.Seq import Seq


def main():
    """ runs main script """

    # Takes standard input and converts to Seq Object to generate
    # reverse complement
    std_input = [Seq(x.rstrip()) for x in sys.stdin]

    # forward strand kmer + 1
    forward_strand = [str(x) for x in std_input]
    # reverse strand kmer + 1
    reverse_strand = [str(r.reverse_complement()) for r in std_input]

    # define all nodes as the set union of forward and reverse strands
    nodes = set(forward_strand) | set(reverse_strand)

    # edges defined as the left and right kmers
    edges = [(kmer[:-1], kmer[1:]) for kmer in nodes]

    with open('deb-output.txt', 'w') as output:
        for edge_pair in sorted(edges):
            output.write("({}, {})\n".format(edge_pair[0], edge_pair[1]))

if __name__ == '__main__':
    main()
