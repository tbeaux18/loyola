#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 01-22-2019

set_operations.py

"""

import sys


def convert2set(std_input_str):
    """ takes std_input string and converts to a real set of numbers """
    std_input_str = std_input_str.replace('{', '')
    std_input_str = std_input_str.replace('}', '')
    new_set = set([int(s) for s in std_input_str.split(',')])

    return new_set

def main():
    """ runs main script """

    std_input = [i.rstrip() for i in sys.stdin]

    pos_int = int(std_input[0])

    set_A = convert2set(std_input[1])

    set_B = convert2set(std_input[2])

    set_union = set_A | set_B

    set_intersection = set_A & set_B

    set_diff_A = set_A - set_B

    set_diff_B = set_B - set_A

    universe = set([i for i in range(1, pos_int)])

    set_u_A = universe - set_A

    set_u_B = universe - set_B


    with open('set-output.txt', 'w') as output:
        output.write(str(set_union)+'\n'+\
        str(set_intersection)+'\n'+\
        str(set_diff_A)+'\n'+\
        str(set_diff_B)+'\n'+\
        str(set_u_A)+'\n'+\
        str(set_u_B))

if __name__ == '__main__':
    main()
