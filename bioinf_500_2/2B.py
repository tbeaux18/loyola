#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-11-2019

2B.py

B) Construct the overlap graph of a collection of k-mers.
Input: A collection Patterns of k-mers.
Output: The overlap graph OVERLAP(Patterns)
The overlap graph can be represented as either an adjacency matrix or an adjacency list.

Input:
    kmer pattern 1 \n
    kmer pattern n \n
Output:
    Kmer is adjact to [List of adjacent sequences]
"""

import sys


class Node:
    """ Node object that contains record id and sequence information
        from the FASTA file. Kmer is set to 3, but can be changedself.
        Attributes:
            record_id (str)
            sequence (str)
            kmer (int)
            connected_nodes (lst)
        Methods:
            add_connected_nodes(self, conn_node)
            get_conn_nodes(self)
            get_record_id(self)
            get_node_sequence(self)
            get_suffix(self)
            get_prefix(self)
    """

    def __init__(self, sequence):
        self.node = sequence
        self.kmer = int(1)
        self.connected_nodes = []

    def add_connected_nodes(self, conn_node):
        """ adds the connected record ID to connected_nodes """
        self.connected_nodes.append(conn_node)

    def get_conn_nodes(self):
        """ returns connected nodes """

        return self.connected_nodes

    def get_node_sequence(self):
        """ returns sequence """

        return self.node

    def get_suffix(self):
        """ returns kmer length suffix of sequence """

        return self.node[-(self.kmer):]

    def get_prefix(self):
        """ returns kmer length prefix of sequence """

        return self.node[:self.kmer]


class OverlapGraph:
    """ Overlap Graph Object that consists of many Node Objects
        Attributes:
            node_dict (dct)
            node_prefix (dct)
            num_of_nodes (int)
        Methods:
            add_node(self, record_id, sequence)
            add_edges(self)
            get_adj_list(self, rec_id)
     """

    def __init__(self):
        self.node_dict = {}
        self.node_prefix = {}
        self.num_of_nodes = 0

    def add_node(self, sequence):
        """ Adds a node object to node_dict instance """
        self.num_of_nodes += 1
        new_node = Node(sequence)
        self.node_dict[sequence] = new_node
        self.node_prefix[sequence] = new_node.get_prefix()
        return new_node

    def add_edges(self):
        """ creates the edges in a list format that the suffix matches the prefix """
        for node_value in self.node_dict.values():
            for prefix_key, prefix_value in self.node_prefix.items():
                if node_value.get_suffix() == prefix_value \
                and node_value.get_node_sequence() != prefix_key:
                    node_value.add_connected_nodes(prefix_key)

    def get_adj_list(self, rec_id):
        """ allows for access to the list of connected nodes on the Node Object """
        return self.node_dict[rec_id].get_conn_nodes()



def main():
    """ runs main script """

    # takes in argument, include the .txt or .fa
    kmer_list = sys.stdin.read().splitlines()

    # Initialize an OverlapGraph object
    overlap_graph = OverlapGraph()

    # Creating Node Objects within the Overlap Graph object
    for kmer in kmer_list:
        overlap_graph.add_node(kmer)

    # Adds the edges to the graph connecting the suffix and prefix of each sequence
    # Utilizes 3mers only for the moment
    overlap_graph.add_edges()

    with open('2B-output.txt', 'w') as output:
        for kmer in kmer_list:
            if overlap_graph.get_adj_list(kmer):
                output.write("{} is adjacent to {}.\n".format(\
                kmer, ', '.join(overlap_graph.get_adj_list(kmer))))
            else:
                output.write("{} is adjacent to nothing.\n".format(kmer))

if __name__ == '__main__':
    main()
