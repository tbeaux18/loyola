#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 10/15/2018

"""

PROBLEM_FILE = input("Enter --> ")

def open_read_write(file_to_open=PROBLEM_FILE):
    """ open a text file and return every even line; indexed at 0 """

    with open(file_to_open, 'r') as working_file:
        # turns the opened file into a list using a list comprehension
        working_file = [sentence.strip() for sentence in working_file.readlines()]

        # iterates through the enumerated file displaying a tuple
        # of (index, string)
        for sentence in enumerate(working_file): #iterates through the file
            if sentence[0] % 2: #Checks to see if the index is divisible by 2
                print(sentence[1]) #prints the string


def main():
    """ prints every even line """

    open_read_write()


if __name__ == '__main__':
    main()
