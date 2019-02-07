
dna = ['ACGCGCGTACGCGCGT', 'ACGTTTAACGACGTTTAACG']

for seq in dna:
    seq_gc = (seq.count('C') + seq.count('G')) / len(seq)
    if seq_gc > 0.55:
        print(seq, seq_gc)
    #print(seq, seq.count('C'))
    #print(seq, seq.count('G'))
