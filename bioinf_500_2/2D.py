#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-11-2019

2D.py

D) Construct the de Bruijn graph from a set of k-mers.
Input: A collection of k-mers Patterns.
Output: The de Bruijn graph DEBRUIJN(Patterns)
The de Bruijn graph can be represented as either an adjacency matrix or an adjacency list.

Requires Biopython

Standard input:
    pattern 1 \n
    pattern n \n
Output:
    Kmer: Kmer L-Kmer: kmer[:-1] R-Kmer: kmer[1:]
"""

import sys
import collections
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
    edges = set(forward_strand) | set(reverse_strand)

    # nodes are the key, edges are the values
    # the left and right kmers are nodes
    # the full kmer is the edge
    # if two nodes have the same edge, they are connected
    deb_graph = collections.defaultdict(list)

    for kmer in edges:
        # left kmer
        deb_graph[kmer[:-1]].append(kmer)
        # right kmer
        deb_graph[kmer[1:]].append(kmer)

    with open('2D-output.txt', 'w') as output:
        for key, value in deb_graph.items():
            output.write("Node: {} Edges: {}\n".format(key, value))
if __name__ == '__main__':
    main()
