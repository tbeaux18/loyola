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

    # forward strand kmer + 1
    forward_strand = set([Seq(x.rstrip()) for x in sys.stdin])

    # reverse strand kmer + 1
    reverse_strand = [r.reverse_complement() for r in forward_strand]

    print(reverse_strand)
if __name__ == '__main__':
    main()
