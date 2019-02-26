#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-06-2019

overlap.py

"""

import sys
from Bio import SeqIO


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

    def __init__(self, record_id, sequence):
        self.record_id = record_id
        self.node = sequence
        self.kmer = int(3)
        self.connected_nodes = []

    def add_connected_nodes(self, conn_node):
        """ adds the connected record ID to connected_nodes """
        self.connected_nodes.append(conn_node)

    def get_conn_nodes(self):
        """ returns connected nodes """

        return self.connected_nodes

    def get_record_id(self):
        """ returns record_id """

        return self.record_id

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

    def add_node(self, record_id, sequence):
        """ Adds a node object to node_dict instance """
        self.num_of_nodes += 1
        new_node = Node(record_id, sequence)
        self.node_dict[record_id] = new_node
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
    fasta_text = sys.argv[1]

    # creates a list of SeqIO objects for easy FASTA parsing
    fasta_list = list(SeqIO.parse(fasta_text, "fasta"))


    # Creating the initial dictionary that holds id and sequence
    seq_dict = {}
    for record in fasta_list:
        seq_dict[str(record.seq)] = str(record.id)


    # Initialize an OverlapGraph object
    overlap_graph = OverlapGraph()


    # Creating Node Objects within the Overlap Graph object
    for key, value in seq_dict.items():
        overlap_graph.add_node(value, key)

    # Adds the edges to the graph connecting the suffix and prefix of each sequence
    # Utilizes 3mers only for the moment
    overlap_graph.add_edges()

    # initializes a temp formatting dict to how rosalind wants the result
    format_dict = {}

    for key, value in seq_dict.items():
        result = overlap_graph.get_adj_list(value)
        if result:
            format_dict[value] = result

    with open('new-output.txt', 'w') as output:
        for key, value in format_dict.items():
            for val in value:
                output.write("{} {}\n".format(key, seq_dict[val]))


if __name__ == '__main__':
    main()
