#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 01-27-2019

1B.py

B) Find all occurrences of a pattern in a string.
Input: Strings Pattern and Genome.
Output: All starting positions in Genome where Pattern appears as a substring. (Note, overlapping
occurrences of Pattern in Genome should be reported; both the original and complementary strands
must be examined. Positions on the complementary strand should be denoted as a negative position.)

Dependencies:
    biopython
Input:
    Text file received as standard input (max 3 lines)
        Line 1) DNA pattern to match
        Line 1) DNA sequence to search
Output:
    Text file showing Sense = [positions] & Antisense = [-(positions)]
"""

import sys
import re
import collections
import time
from Bio.Seq import Seq


def find_pattern_position(input_array):
    """ takes standard input list of 1) dna sequence 2) pattern to match
        and returns a dict of all overlapping positions (1-based index) for both sense
        and antisense strands
        Args:
            stndrd_input (list) : len(2) 1) dna_sequence 2) pattern to match
        Returns:
            position_dict (dict) : sense and antisense positions
    """

    # Creating Biopython Seq object to generate reverse complement
    pos_seq = Seq(input_array[1])
    neg_seq = pos_seq.reverse_complement()

    # generating a match object for each sequence
    pos_matchobj = re.finditer('(?={0})'.format(re.escape(input_array[0])), str(pos_seq))
    neg_matchobj = re.finditer('(?={0})'.format(re.escape(input_array[0])), str(neg_seq))

    # temp dict to generate list of start positions
    temp_dict = collections.defaultdict(list)

    # iterates through each match object
    # TODO: No way to boolean check regex match object.
    # Generates an iterator regardless of True or False
    # Temporary fix in place
    for pos_position in pos_matchobj:
        temp_dict['Sense'].append(pos_position.start()+1)

    for neg_position in neg_matchobj:
        temp_dict['Antisense'].append((neg_position.start()+1)*-1)

    # initializing final dict to return
    position_dict = {}

    # hacked way of checking whether the regex match object contains positions
    for key, value in temp_dict.items():
        if len(temp_dict) < 2:
            if key == 'Sense':
                position_dict[key] = value
                position_dict['Antisense'] = "Pattern not found in antisense sequence."
            else:
                position_dict[key] = value
                position_dict['Sense'] = "Pattern not found in sense sequence."
        else:
            position_dict[key] = value

    return position_dict



def main():
    """ takes in standard input, if standard input is of appropriate length,
        it writes and runs the find pattern function, otherwise will print to console
        to fix input formatting.
    """

    formatted_array = []

    lines = sys.stdin.read().splitlines()

    if len(lines) > 2:
        formatted_array.append(lines[0])
        formatted_array.append(''.join(lines[1:]))

    elif len(lines) == 2:
        formatted_array.append(lines[0])
        formatted_array.append(lines[1])

    with open('1B-output.txt', 'w') as output:
        for key, value in find_pattern_position(formatted_array).items():
            output.write('{} = {}\n'.format(key, value))

if __name__ == '__main__':
    START_TIME = time.time()
    main()
    print(("--- %s seconds ---" % (time.time() - START_TIME)))
