#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 03-02-2019

readfilter_quality.py

"""


import sys
from Bio import SeqIO



def main():
    """ run main script """

    # std input and retains the newline characters
    std_input = sys.stdin.read().splitlines(True)
    # sets the threshold
    quality_numbers = std_input[0].split()

    threshold = int(quality_numbers[0])
    percentage = int(quality_numbers[1]) / 100


    # writes the rest of the fastq reads to a new file for use in biopython
    with open("fastq_read_text.txt", 'w') as tmp_fastq:
        tmp_fastq.writelines(std_input[1:])

    percent_count = 0

    for rec in SeqIO.parse("fastq_read_text.txt", "fastq"):

        base_count = 0

        for base in rec.letter_annotations["phred_quality"]:

            if base >= threshold:

                base_count += 1

        if float(base_count) / len(rec.letter_annotations["phred_quality"]) >= percentage:

            percent_count += 1

    print(percent_count)



if __name__ == '__main__':
    main()
