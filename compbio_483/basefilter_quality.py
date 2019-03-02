#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 03-02-2019

basefilter_quality.py

"""

import sys
import subprocess


def main():
    # std input and retains the newline characters
    std_input = sys.stdin.read().splitlines(True)

    # sets the threshold
    threshold = int(std_input[0].strip())

    # writes the rest of the fastq reads to a new file for use in biopython
    with open("base_filter_quality.txt", 'w') as tmp_fastq:
        tmp_fastq.writelines(std_input[1:])

    trim_cmd = "java -jar /Users/tim.baker/Downloads/Trimmomatic-0.38/trimmomatic-0.38.jar SE -phred33 {} {} LEADING:{} TRAILING:{}".format("base_filter_quality.txt", "base_output.txt", threshold, threshold)

    subprocess.run(trim_cmd, shell=True)

if __name__ == '__main__':
    main()
