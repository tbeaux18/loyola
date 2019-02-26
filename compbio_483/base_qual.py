#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-24-2019

base_qual.py
"""

import sys
from Bio import SeqIO


def main():
    """ runs main script """

    # std input and retains the newline characters
    std_input = sys.stdin.read().splitlines(True)

    # sets the threshold
    threshold = int(std_input[0].strip())

    # writes the rest of the fastq reads to a new file for use in biopython
    with open("fastq_text.txt", 'w') as tmp_fastq:
        tmp_fastq.writelines(std_input[1:])

    temp_list = []

    for rec in SeqIO.parse("fastq_text.txt", "fastq"):
        temp_list.append(rec.letter_annotations["phred_quality"])

    for i in range(len(temp_list)):
        for i in range(temp)
    # holds the bad reads that do not meet the mean threshold
    # bad_reads = (rec for rec in SeqIO.parse("fastq_text.txt", "fastq") \
    #             for num in rec.letter_annotations["phred_quality"] \
    #             if int(num) < threshold)
    # # counts the reads
    # count = SeqIO.write(bad_reads, "bad_qual.fastq", "fastq")
    #
    # # prints to console
    # print(count)

if __name__ == '__main__':
    main()
