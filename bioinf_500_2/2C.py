#!/usr/bin/env python3
"""
@author: Timothy Baker
@date : 02-11-2019

2C.py

C) Construct the de Bruijn graph of a string.
Input: A string Text and an integer k.
Output: The de Bruijn graph DEBRUIJNk(Text)
The de Bruijn graph can be represented as either an adjacency matrix or an adjacency list.

Input:
    kmer_length (int) \n
    sequence (str) \n
Output:
    Kmer: Kmer L-Kmer: kmer[:-1] R-Kmer: kmer[1:]
"""

import sys
import collections
from Bio.Seq import Seq


def composition_k(sequence, kmer_length):
    """ generates list of kmers of length k and sorts them
        lexographically
        Args:
            sequence (str)
            kmer_lenth (int)
        Returns:
            kmer_list (lst) : sorted lexographically
    """
    kmer_list = []

    for i in range(len(sequence) - kmer_length + 1):
        kmer_list.append(sequence[i:i+kmer_length])

    return sorted(kmer_list)


def main():
    """ runs main script """

    # Takes standard input and converts to Seq Object to generate
    # reverse complement
    std_input = sys.stdin.read().splitlines()

    kmer_length = int(std_input[0])
    sequence = std_input[1]

    if len(std_input) > 2:
        sequence = ''.join(std_input[1:])

    # generation ordered kmers
    kmers_list = composition_k(sequence, kmer_length)

    # convert to seq objects for reverse complement generation
    seqobj_list = [Seq(x) for x in kmers_list]

    # forward strand kmer + 1
    forward_strand = [str(x) for x in seqobj_list]
    # reverse strand kmer + 1
    reverse_strand = [str(r.reverse_complement()) for r in seqobj_list]

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
