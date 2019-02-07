import sys
from Bio import pairwise2
from Bio import SeqIO
from Bio.pairwise2 import format_alignment
filepath = str(sys.argv[1]).strip() #read the input file

#takes the matrix and the sequences s and t as parameters because we are going to
# use those to compute the scoring
def get_alignment(matrix, s, t):
    #set the gap and substitution penalties to their corresponding values
    gap_penalty = 1
    sub_penalty = 1
    s_prime = ""
    t_prime = ""
    i, j = len(s), len(t)
    # we use i and j as positions in the matrix. we get their lengths because we
    # previously initalized the matrix to be size of their lengths. furthermore,
    # we are starting from the position where we had the "score" of the matrix
    while i > 0 and j > 0:
        # if the scores in the matrix are equal and there was a match a the corresponding positions
        # in s and t, then we append that corresponding character to the aligned
        # s and t sequences
        if matrix[i][j] == matrix[i-1][j-1] and s[i-1] == t[j-1]:
            s_prime += s[i-1]
            t_prime += t[j-1]
            i -= 1
            j -= 1
        # look at the spot we are currenty at, but also the spot that is in the same
        # column, but just one row above this looks like this:
        # where we are looking in the same column just one row above
        #          . x
        #          . o
        #
        elif matrix[i][j] == matrix[i-1][j] + gap_penalty:
            s_prime += s[i-1]
            t_prime += '-'
            i -= 1
        # this is the same format as above but we are just looking one column to
        # the left instead of one row above
        elif matrix[i][j] == matrix[i][j-1] + gap_penalty:
            s_prime += '-'
            t_prime += t[j-1]
            j -= 1
        # if s[i-1] != t[j-1]
        # this is the case where otherwise the two sequences do not match at the
        # respective character positions. since this is the case, we move over
        # one and up one (diagonal)
        else:
            s_prime += s[i-1]
            t_prime += t[j-1]
            i -= 1
            j -= 1
    # we want to end at position 1,1, so we get there by just iterating and decrementing
    # until we get there. also, these positions are not a match because they are
    # left over, possibly because the sequences are not the same length
    while i > 0:
        s_prime += s[i-1]
        t_prime += '-'
        i -= 1
    while j > 0:
        s_prime += '-'
        t_prime += t[j-1]
        j -= 1

    # returns the aligned sequences in reverse, because we added them starting
    # from the opposite ends
    return s_prime[::-1], t_prime[::-1]
#using biopython to calculate the minimum edit distance and the alignment.
# what this function does is return a list of tuples of all the possible edits to make
# the two strings the same. Because we want the minimum number of edits, we
# will take the last element of the alignments list because they are ordered from
# most edits to least amount of edits.
sequences = list()
for record in SeqIO.parse(filepath, 'fasta'): #get the sequences from the file
    sequences.append(str(record.seq))

gap_penalty = 1
sub_penalty = 1
s = sequences[0]
t = sequences[1]

s_len = len(s)
t_len = len(t)

#create the 2d array using list comprehension, setting each position to None (empty)
# s will correspond to the columns, t will correspond to the rows
matrix = [[None for i in range(t_len+1)] for j in range(s_len+1)]

#set the first row and the first column to the corresponding index value of the sequence
for row in range(s_len + 1):
    matrix[row][0] = row
for col in range(t_len + 1):
    matrix[0][col] = col

#go through each value in the 2d matrix and set the current value of the position
#of the matrix based off the previous row and column of the matrix. basically
#take the minimum value of the previous element in the row or column and the gap penalty,
# or check and see if the two characters are equal and give it a value of 0 if they
# are the same character
for i in range(1, s_len+1):
        for j in range(1, t_len+1):
            matrix[i][j] = min(matrix[i][j-1] + gap_penalty,
                          matrix[i-1][j] + gap_penalty,
                          matrix[i-1][j-1] + (sub_penalty if s[i-1] != t[j-1] else 0))
#the score will be the last value in the matrix
score = matrix[-1][-1]
s_prime, t_prime = get_alignment(matrix, s, t)
print(score)
print(s_prime)
print(t_prime)
