#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 03-01-2019

assembly_n5075.py

"""
import sys
from itertools import accumulate



def main():
    """ runs main script """

    std_input = sys.stdin.read().splitlines()

    seq_length = sorted([len(x) for x in std_input], reverse=True)
    seq_cumsum = [i for i in accumulate(seq_length)]

    n50_score = int(sum(seq_length)/2)
    n75_score = int(sum(seq_length)*0.75)

    n50_max = max([seq_length[i] for i in range(len(seq_cumsum)) if seq_cumsum[i] >= n50_score])
    n75_max = max([seq_length[i] for i in range(len(seq_cumsum)) if seq_cumsum[i] >= n75_score])

    print(n50_max, n75_max)




if __name__ == '__main__':
    main()
