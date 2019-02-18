#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-11-2019

2E.py

Input:
    kmer pattern 1 \n
    kmer pattern n \n
Output:
    cycle
"""

import sys
import random
from Bio.Seq import Seq


def eulerian_cycle(deb_graph):
    """ no idea if this is correct """
    cycle = []

    random_start_index = random.randint(0, len(deb_graph) - 1)
    random_start_node = deb_graph[random_start_index][0]
    cycle.append(deb_graph[random_start_index][0])

    while len(deb_graph) > 0:
        for edge in deb_graph:
            if edge[0] == random_start_node:
                random_start_node = edge[1]
                deb_graph.remove(edge)
                cycle.append(random_start_node)
                break
            else:
                for edge in deb_graph:
                    random_start_node = edge[0]
                    break
                else:
                    print("no path")

    return cycle


def main():
    """ runs main script """

    # Takes standard input and converts to Seq Object to generate
    # reverse complement
    std_input = [Seq(x.rstrip()) for x in sys.stdin]

    # forward strand kmer + 1
    forward_strand = [str(x) for x in std_input]
    # reverse strand kmer + 1
    reverse_strand = [str(r.reverse_complement()) for r in std_input]

    # define all edges as the set union of forward and reverse strand kmers
    edges = set(forward_strand) | set(reverse_strand)

    graph = [(kmer[:-1], kmer[1:]) for kmer in edges]

    cycle = eulerian_cycle(graph)

    with open('2E-output.txt', 'w') as output:
        output.write(' > '.join(cycle) + '\n')

if __name__ == '__main__':
    main()
