#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 10/10/2018

Dictionary Rosalind Problem 6
"""

import sys


def create_word_dict(sentence):
    """
    Args:
        sentence (str) : stdin string with spaces; required for word_list
    Returns:
        dictionary count
    """
    word_list = sentence.split(' ') # turns string into a list of words

    # key = word, value = count
    word_count_dict = {} #initializes the dictionary

    for word in word_list: #iterates through the list
        #sets word as key and count as value
        word_count_dict[word] = word_list.count(word)

    for key, value in word_count_dict.items(): # iterates through key, value pair
        print('{} {}'.format(key, value)) #pretty prints the key, value




def main():
    """ prints key, value to console """
    line_sentence = sys.stdin.read().rstrip()
    create_word_dict(line_sentence)


if __name__ == '__main__':
    main()
