#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 10/10/2018

String and List Rosalind Problem 3

"""
import sys


def create_word(sentence, a_index, b_index, c_index, d_index):
    """
    Args:
        sentence (str): concatentated sentence
        a_index (int): index 0 of stdin list
        b_index (int): index 1 of stdin list
        c_index (int): index 2 of stdin list
        d_index (int): index 3 of stdin list
    Returns:
        index spliced string(s)
    """
    # Turns string into a list
    sentence_list = list(sentence)

    # splices based off input indices and joins the characters
    first_word = "".join(sentence_list[a_index:b_index+1])
    second_word = "".join(sentence_list[c_index:d_index+1])

    # passes the newly joined string with an inclusive whitespace
    return '%s %s' % (first_word, second_word)


def main():
    """ Reads standard input and prints create_word function to console """

    lines = [line.rstrip() for line in sys.stdin.readlines()]
    experiment_sentence = str(lines[0])
    indices = [int(x) for x in lines[1].split()]
    print(create_word(experiment_sentence, indices[0], indices[1], indices[2], indices[3]))


if __name__ == '__main__':
    main()
