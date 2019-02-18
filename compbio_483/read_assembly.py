#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-16-2019

read_assembly.py

Dependencies:
    networkx
    biopython
"""


import sys
import collections
from random import randint
from Bio.Seq import Seq
import networkx as nx

def generate_iterating_kmer(dna_sequence_list, kmer_length):
    """ Takes a DNA sequence and generates all kmer substrings of the sequence
        Args:
            dna_sequence (str) : must be ACGT; no RNA transcripts
            kmer_length (int) : default is set to 3 for codon, but can be changed
        Returns:
            a list of all possible kmers
    """

    kmer_list = [] # initializes the kmer list; order is important

    for seq in dna_sequence_list:
        for i in range(len(seq) - kmer_length + 1):
            kmer_list.append(seq[i:i+kmer_length])

    return kmer_list


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
    all_reads = list(set(forward_strand) | set(reverse_strand))

    # counts the total number of reads for both forward and reverse
    num_of_reads = len(forward_strand)

    # starts with longest k-length iterator and begins to decrease so long as
    # no solution exists
    for k_length in range(len(all_reads[0]), 2, -1):

        # creates the kmers of k-length from all reads
        all_reads = set(generate_iterating_kmer(all_reads, k_length))

        # initializes a de bruijnized graph
        deb_graph = collections.defaultdict(list)

        # creates the de bruijn graph for the left and right k-1-mers
        for kmer in all_reads:
            deb_graph[kmer[:-1]].append(kmer[1:])
            deb_graph[kmer[1:]].append(kmer[:-1])

        # initializes an instance of Directed Graph from NetworkX python package
        # deb_graph provides the node and edge coordinates
        directed_graph = nx.DiGraph(deb_graph)

        # initializes a copy of the graph along with all nodes and edges
        copy_digraph = directed_graph.copy()

        # iterates through each node on the graph, and if the node has a
        # degree of 2 (signalling 2 directed cycles and even euler path)
        # adds those nodes to the copy_digraph to cycle through
        for a_node in directed_graph.__iter__():
            if directed_graph.degree(a_node) == 2:
                copy_digraph.add_node(a_node)

        # begins to find all the cycles of the copy di_graph, since all nodes
        # have 2 degrees
        cycle_list = list(nx.simple_cycles(copy_digraph))

        # initialize an empty list to store superstrings
        list_of_superstrings = []

        # iterates through each of the cycles, and only considers cycles that
        # have equal or greater amount of elements per cycle
        # this will indicate that each node is in each of the starting
        # reads at least once. Also takes the first letter of each element of the
        # cycle list to construct the sequence
        for cycle in cycle_list:
            if len(cycle) >= num_of_reads:
                list_of_superstrings.append(''.join(x[0] for x in cycle if x in ''.join(all_reads)))

        # list will be empty until it iterates to a kmer length that can be cyclic
        # randomly choosing a sequence from the list to output
        # break must be included otherwise it will iterate all the way down to
        # 2 kmers which could present incorrect assemblies
        if list_of_superstrings:
            i = randint(0, len(list_of_superstrings) - 1)
            with open('assembly-output.txt', 'w') as output:
                output.write(list_of_superstrings[i])
            break



if __name__ == '__main__':
    main()
