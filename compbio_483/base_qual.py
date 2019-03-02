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

    record_list = []
    for rec in SeqIO.parse("fastq_text.txt", "fastq"):
        record_list.append(rec.letter_annotations["phred_quality"])

    per_base_quality = [float(0) for x in range(len(record_list[0]))]

    for record in record_list:
        for idx in enumerate(record):
            per_base_quality[idx[0]] += record[idx[0]]

    for index in enumerate(per_base_quality):
        per_base_quality[index[0]] = per_base_quality[index[0]] / len(record_list)

    base_numbers = len([i for i in per_base_quality if i < threshold])

    print(base_numbers)


if __name__ == '__main__':
    main()
