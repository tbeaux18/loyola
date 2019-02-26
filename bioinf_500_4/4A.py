#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-25-2019

4A.py

Input:
    flu_fasta.fa
    rhino_fasta.fa
    patient_fastq.txt

When running the script, order is important:
    1) flu fasta 2) rhino fasta 3) patient fastq
    python3 4A.py influenza_ref.fasta rhinovir_ref.fasta patient_fastq.txt

"""

import sys
from suffix_trees import STree
from Bio import SeqIO


def main():
    """ runs main script """

    input_files = sys.argv

    flu_fasta = input_files[1]
    rhino_fasta = input_files[2]
    patient_fastq = input_files[3]

    flu_fasta_ref = list(SeqIO.parse(flu_fasta, "fasta"))
    rhino_fasta_ref = list(SeqIO.parse(rhino_fasta, "fasta"))
    patient_fastq_reads = list(SeqIO.parse(patient_fastq, "fastq"))

    flu_seq = [str(flu.seq) for flu in flu_fasta_ref]
    rhi_seq = [str(rhi.seq) for rhi in rhino_fasta_ref]

    flu_suffix = STree.STree(flu_seq)
    rhi_suffix = STree.STree(rhi_seq)

    flu_score = 0
    rhi_score = 0

    positive_strain_list = []

    for read in patient_fastq_reads:

        if flu_suffix.find_all(str(read.seq)):
            flu_score += 1
            for record in flu_fasta_ref:
                if str(read.seq) in record.seq:
                    positive_strain_list.append(str(record.description))

        elif rhi_suffix.find_all(str(read.seq)):
            rhi_score += 1
            for record in rhino_fasta_ref:
                if str(read.seq) in record.seq:
                    positive_strain_list.append(str(record.description))


    with open('patient-report.txt', 'w') as output:
        output.write("Sequencing results are detecting the following:\n")
        output.write('\n'.join(positive_strain_list))
        if flu_score and not rhi_score:
            output.write("\nPatient is positive for the Influenza A virus.")
        elif rhi_score and not flu_score:
            output.write("\nPatient is positive for the Human Rhinovirus Strain 89")
        elif not flu_score and not rhi_score:
            output.write("\nInfluenza A virus and Human Rhinovirus not detected in patient.")



if __name__ == '__main__':
    main()
