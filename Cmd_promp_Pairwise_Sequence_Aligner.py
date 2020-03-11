
"""
Created on Fri Jan 31 11:16:25 2020
Final edits on March 10 2020

@author: crlsalmeida
"""
# Analysis Preparation========================================================================================


print("\n")
print("\n")

print("Greeting.")
print("This is a simple Pairwise Sequence Aligner created by crlsalmeida.")
print("It is a simple aligner because there are no affine gap penalties in this aligner.")
print("\n")
print("Before getting started, move the 2 FASTA files you wish to align to the same folder as this Aligner.")
print("Let's get started.")

print("\n")
print("\n")

# Import libraries

import numpy as np
import pandas as pd
import os

# Run the following function reads and retrive the sequences from fasta files as strings

def read_fasta(fasta1, fasta2):
    
    '''
    This function opens 2 .fasta files
    It will return the first line of both .fasta files as 'title1' and 'title2'
    It will return both sequences as strings
    
    '''
    # os.chdir(workingdirectory)
    s1 = ''
    s2 = ''
    title1 = ''
    title2 = ''
    with open(fasta1, 'r') as FASTA1:
        for line in FASTA1.readlines():
            if line[0] == '>':
                title1 = line
            elif not line[0] == '>':
                s1 += line.rstrip()
    with open(fasta2, 'r') as FASTA2:
        for line in FASTA2.readlines():
            if line[0] == '>':
                title2 = line
            elif not line[0] == '>':
                s2 += line.rstrip()
    return s1,s2, title1, title2

# Command line call

fasta1 = input("Please enter the file name of the first sequence (ex. sequence1.fa): ")
print("\n")
fasta2 = input("Please enter the file name of the second sequence (ex. sequence2.fa): ")

s1 = read_fasta(fasta1, fasta2)[0]
s2 = read_fasta(fasta1, fasta2)[1]
titleseq1 = read_fasta(fasta1, fasta2)[2]
titleseq2 = read_fasta(fasta1, fasta2)[3]

# Score Matrices and Look-up Functions==================================================================================

# The following functions will be used in the Needleman-Wunsch and Smith-Waterman Functions
    
# Below is the functioned that corresponds to aligning nucleotide sequences
    
def NucSeq (ch1, ch2, ScrMtrx):
    
    '''
    
    This is a look-up function that corresponds to Match and Mismatch scores
    It will return the the row and column location of nucleotides
        
    '''
    
    if ch1 == 'A':
        idx1 = 0
    elif ch1 == 'C':
        idx1 = 1
    elif ch1 == 'G':
        idx1 = 2
    else:
        idx1 = 3
    if ch2 == 'A':
        idx2 = 0
    elif ch2 == 'C':
        idx2 = 1
    elif ch2 == 'G':
        idx2 = 2
    else:
        idx2 = 3
    return idx1, idx2

# Below is the function that corresponds to alignment of amino acid sequences

def AASeq (ch1, ch2, ScrMtrx):
    
    '''
    
    This is a look-up function that corresponds PAM or BLOSUM tables
    
    In order to function properly, Amino Acids have to be entered into the score tables in the following order:
    A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  * 
    
    '''
    
    if ch1 == 'A':
        idx1 = 0
    elif ch1 == 'R':
        idx1 = 1
    elif ch1 == 'N':
        idx1 = 2
    elif ch1 == 'D':
        idx1 = 3
    elif ch1 == 'C':
        idx1 = 4
    elif ch1 == 'Q':
        idx1 = 5
    elif ch1 == 'E':
        idx1 = 6
    elif ch1 == 'G':
        idx1 = 7
    elif ch1 == 'H':
        idx1 = 8
    elif ch1 == 'I':
        idx1 = 9
    elif ch1 == 'L':
        idx1 = 10
    elif ch1 == 'K':
        idx1 = 11
    elif ch1 == 'M':
        idx1 = 12
    elif ch1 == 'F':
        idx1 = 13
    elif ch1 == 'P':
        idx1 = 14
    elif ch1 == 'S':
        idx1 = 15
    elif ch1 == 'T':
        idx1 = 16
    elif ch1 == 'W':
        idx1 = 17
    elif ch1 == 'Y':
        idx1 = 18
    elif ch1 == 'V':
        idx1 = 19
    elif ch1 == 'B':
        idx1 = 20
    elif ch1 == 'Z':
        idx1 = 21
    elif ch1 == 'X':
        idx1 = 22
    elif ch1 == '*':
        idx1 = 23
    if ch2 == 'A':
        idx2 = 0
    elif ch2 == 'R':
        idx2 = 1
    elif ch2 == 'N':
        idx2 = 2
    elif ch2 == 'D':
        idx2 = 3
    elif ch2 == 'C':
        idx2 = 4
    elif ch2 == 'Q':
        idx2 = 5
    elif ch2 == 'E':
        idx2 = 6
    elif ch2 == 'G':
        idx2 = 7
    elif ch2 == 'H':
        idx2 = 8
    elif ch2 == 'I':
        idx2 = 9
    elif ch2 == 'L':
        idx2 = 10
    elif ch2 == 'K':
        idx2 = 11
    elif ch2 == 'M':
        idx2 = 12
    elif ch2 == 'F':
        idx2 = 13
    elif ch2 == 'P':
        idx2 = 14
    elif ch2 == 'S':
        idx2 = 15
    elif ch2 == 'T':
        idx2 = 16
    elif ch2 == 'W':
        idx2 = 17
    elif ch2 == 'Y':
        idx2 = 18
    elif ch2 == 'V':
        idx2 = 19
    elif ch2 == 'B':
        idx2 = 20
    elif ch2 == 'Z':
        idx2 = 21
    elif ch2 == 'X':
        idx2 = 22
    elif ch2 == '*':
        idx2 = 23
    return idx1, idx2

PAM250 = np.array([[ 2, -2,  0,  0, -2,  0,  0,  1, -1, -1, -2, -1, -1, -3,  1,  1,  1, -6, -3,  0,  0,  0,  0, -8],
                   [-2,  6,  0, -1, -4,  1, -1, -3,  2, -2, -3,  3,  0, -4,  0,  0, -1,  2, -4, -2, -1,  0, -1, -8],
                   [ 0,  0,  2,  2, -4,  1,  1,  0,  2, -2, -3,  1, -2, -3,  0,  1,  0, -4, -2, -2,  2,  1,  0, -8],
                   [ 0, -1,  2,  4, -5,  2,  3,  1,  1, -2, -4,  0, -3, -6, -1,  0,  0, -7, -4, -2,  3,  3, -1, -8],
                   [-2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3,  0, -2, -8,  0, -2, -4, -5, -3, -8],
                   [ 0,  1,  1,  2, -5,  4,  2, -1,  3, -2, -2,  1, -1, -5,  0, -1, -1, -5, -4, -2,  1,  3, -1, -8],
                   [ 0, -1,  1,  3, -5,  2,  4,  0,  1, -2, -3,  0, -2, -5, -1,  0,  0, -7, -4, -2,  3,  3, -1, -8],
                   [ 1, -3,  0,  1, -3, -1,  0,  5, -2, -3, -4, -2, -3, -5,  0,  1,  0, -7, -5, -1,  0,  0, -1, -8],
                   [-1,  2,  2,  1, -3,  3,  1, -2,  6, -2, -2,  0, -2, -2,  0, -1, -1, -3,  0, -2,  1,  2, -1, -8],
                   [-1, -2, -2, -2, -2, -2, -2, -3, -2,  5,  2, -2,  2,  1, -2, -1,  0, -5, -1,  4, -2, -2, -1, -8],
                   [-2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6, -3,  4,  2, -3, -3, -2, -2, -1,  2, -3, -3, -1, -8],
                   [-1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5,  0, -5, -1,  0,  0, -3, -4, -2,  1,  0, -1, -8],
                   [-1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6,  0, -2, -2, -1, -4, -2,  2, -2, -2, -1, -8],
                   [-3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,  9, -5, -3, -3,  0,  7, -1, -4, -5, -2, -8],
                   [ 1,  0,  0, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2, -5,  6,  1,  0, -6, -5, -1, -1,  0, -1, -8],
                   [ 1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2, -3,  1,  2,  1, -2, -3, -1,  0,  0,  0, -8],
                   [ 1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1, -3,  0,  1,  3, -5, -3,  0,  0, -1,  0, -8],
                   [-6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,  0, -6, -2, -5, 17,  0, -6, -5, -6, -4, -8],
                   [-3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,  7, -5, -3, -3,  0, 10, -2, -3, -4, -2, -8],
                   [ 0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2, -1, -1, -1,  0, -6, -2,  4, -2, -2, -1, -8],
                   [ 0, -1,  2,  3, -4,  1,  3,  0,  1, -2, -3,  1, -2, -4, -1,  0,  0, -5, -3, -2,  3,  2, -1, -8],
                   [ 0,  0,  1,  3, -5,  3,  3,  0,  2, -2, -3,  0, -2, -5,  0,  0, -1, -6, -4, -2,  2,  3, -1, -8],
                   [ 0, -1,  0, -1, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1,  0,  0, -4, -2, -1, -1, -1, -1, -8],
                   [-8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1]], dtype = float)

BLOSUM62 = np.array([[ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4], 
                     [-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4], 
                     [-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4], 
                     [-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4], 
                     [ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4], 
                     [-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4], 
                     [-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4], 
                     [ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4], 
                     [-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4], 
                     [-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4], 
                     [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4], 
                     [-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4], 
                     [-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4], 
                     [-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4], 
                     [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4], 
                     [ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4], 
                     [ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4], 
                     [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4], 
                     [-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4], 
                     [ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4], 
                     [-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4], 
                     [-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4], 
                     [ 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4], 
                     [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1]], dtype = float)


BLOSUM80 = np.array([[ 7, -3, -3, -3, -1, -2, -2,  0, -3, -3, -3, -1, -2, -4, -1,  2,  0, -5, -4, -1, -3, -2, -1, -8],
                     [-3,  9, -1, -3, -6,  1, -1, -4,  0, -5, -4,  3, -3, -5, -3, -2, -2, -5, -4, -4, -2,  0, -2, -8],
                     [-3, -1,  9,  2, -5,  0, -1, -1,  1, -6, -6,  0, -4, -6, -4,  1,  0, -7, -4, -5,  5, -1, -2, -8],
                     [-3, -3,  2, 10, -7, -1,  2, -3, -2, -7, -7, -2, -6, -6, -3, -1, -2, -8, -6, -6,  6,  1, -3, -8],
                     [-1, -6, -5, -7, 13, -5, -7, -6, -7, -2, -3, -6, -3, -4, -6, -2, -2, -5, -5, -2, -6, -7, -4, -8],
                     [-2,  1,  0, -1, -5,  9,  3, -4,  1, -5, -4,  2, -1, -5, -3, -1, -1, -4, -3, -4, -1,  5, -2, -8],
                     [-2, -1, -1,  2, -7,  3,  8, -4,  0, -6, -6,  1, -4, -6, -2, -1, -2, -6, -5, -4,  1,  6, -2, -8],
                     [ 0, -4, -1, -3, -6, -4, -4,  9, -4, -7, -7, -3, -5, -6, -5, -1, -3, -6, -6, -6, -2, -4, -3, -8],
                     [-3,  0,  1, -2, -7,  1,  0, -4, 12, -6, -5, -1, -4, -2, -4, -2, -3, -4,  3, -5, -1,  0, -2, -8],
                     [-3, -5, -6, -7, -2, -5, -6, -7, -6,  7,  2, -5,  2, -1, -5, -4, -2, -5, -3,  4, -6, -6, -2, -8],
                     [-3, -4, -6, -7, -3, -4, -6, -7, -5,  2,  6, -4,  3,  0, -5, -4, -3, -4, -2,  1, -7, -5, -2, -8],
                     [-1,  3,  0, -2, -6,  2,  1, -3, -1, -5, -4,  8, -3, -5, -2, -1, -1, -6, -4, -4, -1,  1, -2, -8],
                     [-2, -3, -4, -6, -3, -1, -4, -5, -4,  2,  3, -3,  9,  0, -4, -3, -1, -3, -3,  1, -5, -3, -2, -8],
                     [-4, -5, -6, -6, -4, -5, -6, -6, -2, -1,  0, -5,  0, 10, -6, -4, -4,  0,  4, -2, -6, -6, -3, -8],
                     [-1, -3, -4, -3, -6, -3, -2, -5, -4, -5, -5, -2, -4, -6, 12, -2, -3, -7, -6, -4, -4, -2, -3, -8],
                     [ 2, -2,  1, -1, -2, -1, -1, -1, -2, -4, -4, -1, -3, -4, -2,  7,  2, -6, -3, -3,  0, -1, -1, -8],
                     [ 0, -2,  0, -2, -2, -1, -2, -3, -3, -2, -3, -1, -1, -4, -3,  2,  8, -5, -3,  0, -1, -2, -1, -8],
                     [-5, -5, -7, -8, -5, -4, -6, -6, -4, -5, -4, -6, -3,  0, -7, -6, -5, 16,  3, -5, -8, -5, -5, -8],
                     [-4, -4, -4, -6, -5, -3, -5, -6,  3, -3, -2, -4, -3,  4, -6, -3, -3,  3, 11, -3, -5, -4, -3, -8],
                     [-1, -4, -5, -6, -2, -4, -4, -6, -5,  4,  1, -4,  1, -2, -4, -3,  0, -5, -3,  7, -6, -4, -2, -8],
                     [-3, -2,  5,  6, -6, -1,  1, -2, -1, -6, -7, -1, -5, -6, -4,  0, -1, -8, -5, -6,  6,  0, -3, -8],
                     [-2,  0, -1,  1, -7,  5,  6, -4,  0, -6, -5,  1, -3, -6, -2, -1, -2, -5, -4, -4,  0,  6, -1, -8],
                     [-1, -2, -2, -3, -4, -2, -2, -3, -2, -2, -2, -2, -2, -3, -3, -1, -1, -5, -3, -2, -3, -1, -2, -8],
                     [-8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1]], dtype = float)


# Command line call
print("\n")
typeofsequence = input("Are you aligning (1) DNA or (2) Amino Acid sequences? Enter either 1 or 2: ")

if typeofsequence == '1':
    SeqType = NucSeq
    print("\n")
    print("DNA - understood.")
elif typeofsequence == '2':
    SeqType = AASeq
    print("\n")
    print("Amino Acid - got it.")
    print("Due to my programing limitations, I can only offer you three options -")
else:
    SeqType = NucSeq
    print("\n")
    print("You're pretty bad at following intructions aren't you?")
    print("Anyway, I'm going to assume you meant to say (1) DNA.")
    print("If you don't want to align two DNA sequences, then you have to restart this script from the beginning.")



if SeqType == NucSeq:
    print("\n")
    MatchScore = input("Please specify the Match Score: ")
    print("\n")
    MismatchScore = input("Please specify the Mismatch Score: ")
    print("\n")
    GapScore = input("Please specify the Gap Score: ")
    NucScrTbl = np.array([[MatchScore, MismatchScore, MismatchScore, MismatchScore],
                      [MismatchScore, MatchScore, MismatchScore, MismatchScore],
                      [MismatchScore, MismatchScore, MatchScore, MismatchScore],
                      [MismatchScore, MismatchScore, MismatchScore, MatchScore]], dtype = float)
    ScrMtrx = NucScrTbl
    MatrixName = 'name'
else:
    print("\n")
    AAMatrix = input("Would you like to use a (1) PAM250, (2) BLOSUM62, or (3) BLOSUM80 score matrix? Enter 1, 2, or 3: ")
    if AAMatrix == '1':
        ScrMtrx = PAM250
        MatrixName = "PAM250"
        MatchScore = 0
        MismatchScore = 0
        print("\n")
        print("PAM250 - cool beans.")
    elif AAMatrix == '2':
        SeqType = BLOSUM62
        MatrixName = "BLOSUM62"
        MatchScore = 0
        MismatchScore = 0
        print("\n")
        print("BLOSUM62 - good choice, dude.")
    elif AAMatrix == '3':
        ScrMtrx = BLOSUM80
        MatrixName = "BLOSUM80"
        MatchScore = 0
        MismatchScore = 0
        print("\n")
        print("BLOSUM80 - I like your style.")
    else:
        ScrMtrx = PAM250
        MatrixName = "PAM250"
        MatchScore = 0
        MismatchScore = 0
        print("\n")
        print("You're pretty bad at following intructions aren't you?")
        print("Anyway, I'm going to assume you meant to say (1) PAM250.")
        print("If you don't want to use a PAM250 matrix, then you have to restart this script from the beginning.")
    print("\n")
    GapScore = input("Please specify the Gap Score: ")

# Needleman-Wunsch Function==================================================================================

def NeedlemanWunsch(s1,s2, GapScore, ScrMtrx, SeqType):
    
    '''
    
    This function performs a Needleman-Wunsch global sequence alignment without affine gap penalty
    By default, this is set up to align nucleotide sequences, SeqType
        Define the SeqType as either 'NucSeq' for nucleotide sequences
        Or 'AASeq' for amino acid sequences
    The default scoring matrix, ScrMtrx, is a nucleotide score matrix
        Define ScrMtrx as either 'NucScrTbl' for nucleotide sequences
        Or 'PAM250', 'BLOSUM62', etc. for AA sequences provided the tables were loaded as described in the 'AASeq' function

    '''
    
    # Set up values for Matrices
    
    n = len(s1) # n is associated with x, i
    m = len(s2) # m is associated with y, j
    s1L = list(s1.strip()) # convert s1 string to list
    s2L = list(s2.strip()) # convert s2 string to list
    d = float(GapScore)
    
    # Create Score Matrix
    
    M = np.zeros((m+1, n+1)) # create Score Matrix

    # Assign values for first row and column
    for i in range(0,len(s1)+1):
        M[0][i] = i*d
    for j in range(0, len(s2)+1):
        M[j][0] = j*d
            
    # Create Direction Matrix
    
    D = np.zeros((m+1, n+1))
    D[0,1:] = 2 # denotes vertical/gap
    D[1:,0] = 3 # denotes horizontal/gap   
    
    # Calculate Values of Score Matrix and enter directions in Direction Matrix
        
    for i in range(1,len(s1)+1): # Note that here, x corresponds to s1
        for j in range(1,len(s2)+1): # Note that here, i corresponds to s2
            ch1 = s1[i-1]
            ch2 = s2[j-1]
            rowidx = SeqType(ch1, ch2, ScrMtrx)[0]
            colidx = SeqType(ch1, ch2, ScrMtrx)[1]
            S = ScrMtrx[rowidx, colidx]
            Match = M[j-1][i-1] + S
            Delete = M[j-1][i] + d
            Insert = M[j][i-1] + d
            M[j,i] = max(Match, Delete, Insert)
            if M[j,i] == Match:
                D[j,i] = 1 # denotes diagonal
            elif M[j,i] == Delete:
                D[j,i] = 3 # denotes vertical
            else:
                D[j,i] = 2 # denotes horizontal
                
    # Return alignment
    
    As1 = '' # This is the stand in for alignment of s1
    As2 = '' # This is the stand in for alignment of s2
    Alg = '' # This is the stand in for alignment labels
    Alsp = '' # This is the stand in for space between alignments used for printing, if needed
    
    i = len(s1)
    j = len(s2)
    
    AlignmentScore = M[j,i]
    
    while i>0 or j>0: # Note that j comes before i
        if D[j,i] == 1:
            ch1 = s1[i-1]
            ch2 = s2[j-1]
            rowidx = SeqType(ch1, ch2, ScrMtrx)[0]
            colidx = SeqType(ch1, ch2, ScrMtrx)[1]
            S = ScrMtrx[rowidx, colidx]
            if s1[i-1]==s2[j-1]:
                Alg = "|"+Alg
            elif SeqType is AASeq and S > 0:
                Alg = ":"+Alg
            elif SeqType is AASeq and S <= 0:
                Alg = "."+Alg
            else:
                Alg = " "+Alg
            Alsp = " "+Alsp
            As1 = s1[i-1]+As1
            As2 = s2[j-1]+As2
            i = i-1
            j = j-1
        elif D[j,i] == 2:
            As1 = s1[i-1]+As1
            As2 = "-"+As2
            Alg = " "+Alg
            Alsp = " "+Alsp
            i = i-1
        elif D[j,i] == 3:
            As1 = "-"+As1
            As2 = s2[j-1]+As2
            Alg = " "+Alg
            Alsp = " "+Alsp
            j = j-1
        else:
            break
    
    # Labeling Matrices
        
    s1L.insert(0,"-") # insert - to [0] of s1 List
    s2L.insert(0,"-") # insert - to [0] of s2 List
    M = pd.DataFrame(M, columns=s1L, index=s2L)
    D = pd.DataFrame(D, columns=s1L, index=s2L)
    
    # Returns alingned sequences in a nice matrix format
    
    A1List = list(As1)
    A2List = list(Alg)
    A3List = list(As2)

    AlgM = np.matrix([A1List,A2List,A3List])
            
    width = 75
    i1 = [As1[i:i+width] for i in range(0, len(As1), width)]
    i2 = [Alg[i:i+width] for i in range(0, len(Alg), width)]
    i3 = [As2[i:i+width] for i in range(0, len(As2), width)]
    i4 = [Alsp[i:i+width] for i in range(0, len(Alsp), width)]
    out = [None]*(len(i1)+len(i2)+len(i3)+len(i4))
    out [::4] = i1
    out [1::4] = i2
    out [2::4] = i3
    out [3::4] = i4
        
    # output = ('\n'.join(map(str, out)))
            
    return M,D,As1,Alg,As2,Alsp,AlignmentScore,AlgM,out

# Smith-Waterman Function==================================================================================
    

def SmithWaterman(s1,s2, GapScore, ScrMtrx, SeqType):
    
    '''
    
    This function performs a Smith-Waterman local alignment without affine gap penalty
    By default, this is set up to align nucleotide sequences, SeqType
        Define the SeqType as either 'NucSeq' for nucleotide sequences
        Or 'AASeq' for amino acid sequences
    The default scoring matrix, ScrMtrx, is a nucleotide score matrix
        Define ScrMtrx as either 'NucScrTbl' for nucleotide sequences
        Or 'PAM250', 'BLOSUM62', etc. for AA sequences provided the tables were loaded as described in the 'AASeq' function

    '''
    # Set up values for Matrices
    
    n = len(s1) # n is associated with x, i
    m = len(s2) # m is associated with y, j
    s1L = list(s1.strip()) # convert s1 string to list
    s2L = list(s2.strip()) # convert s2 string to list
    d = float(GapScore)

        
    # Create Score Matrix
    
    M = np.zeros((m+1, n+1)) # create Score Matrix

    # Assign values for first row and column
    for i in range(0,len(s1)+1):
        M[0][i] = i*0
    for j in range(0, len(s2)+1):
        M[j][0] = j*0
            
    # Create Direction Matrix
    
    D = np.zeros((m+1, n+1))
    D[0,1:] = 2 # denotes vertical/gap
    D[1:,0] = 3 # denotes horizontal/gap   
    
    # Calculate Values of Score Matrix and enter directions in Direction Matrix
    
    max_score = 0
    max_i = 0
    max_j = 0
    
    for i in range(1,len(s1)+1):
        for j in range(1,len(s2)+1):
            ch1 = s1[i-1]
            ch2 = s2[j-1]
            rowidx = SeqType(ch1, ch2, ScrMtrx)[0]
            colidx = SeqType(ch1, ch2, ScrMtrx)[1]
            S = ScrMtrx[rowidx, colidx]
            Match = M[j-1][i-1] + S
            Delete = M[j-1][i] + d
            Insert = M[j][i-1] + d
            M[j,i] = max(Match, Delete, Insert,0)
            score = max(Match, Delete, Insert,0)
            if score > max_score:
                max_score = score
                max_i = [i]
                max_j = [j]
            if M[j,i] == Match:
                D[j,i] = 1 # denotes diagonal
            elif M[j,i] == Delete:
                D[j,i] = 3 # denotes vertical
            elif M[j,i] == 0:
                D[j,i] == 0
            else:
                D[j,i] = 2 # denotes horizontal
                
    # Return alignment
    
    As1 = '' # This is the stand in for alignment of s1
    As2 = '' # This is the stand in for alignment of s2
    Alg = '' # This is the stand in for alignment labels
    Alsp = '' # This is the stand in for space between alignments used for printing, if needed
    
    i = max_i[0]
    j = max_j[0]
    
    AlignmentScore = M[j,i]
    
    while i>0 or j>0: # Note that j comes before i
        if M[j,i] <=0:
            break
        elif D[j,i] == 1:
            ch1 = s1[i-1]
            ch2 = s2[j-1]
            rowidx = SeqType(ch1, ch2, ScrMtrx)[0]
            colidx = SeqType(ch1, ch2, ScrMtrx)[1]
            S = ScrMtrx[rowidx, colidx]
            Alsp = " "+Alsp
            As1 = s1[i-1]+As1
            As2 = s2[j-1]+As2
            if s1[i-1]==s2[j-1]:
                Alg = "|"+Alg
            elif SeqType is AASeq and S > 0:
                Alg = ":"+Alg
            elif SeqType is AASeq and S <= 0:
                Alg = "."+Alg
            else:
                Alg = " "+Alg
            i = i-1
            j = j-1
        elif D[j,i] == 2:
            As1 = s1[i-1]+As1
            As2 = "-"+As2
            Alg = " "+Alg
            Alsp = " "+Alsp
            i = i-1
        elif D[j,i] == 3:
            As1 = "-"+As1
            As2 = s2[j-1]+As2
            Alg = " "+Alg
            Alsp = " "+Alsp
            j = j-1
        else:
            break
    
    # Labeling Matrices
        
    s1L.insert(0,"-") # insert - to [0] of s1 List
    s2L.insert(0,"-") # insert - to [0] of s2 List
    M = pd.DataFrame(M, columns=s1L, index=s2L)
    D = pd.DataFrame(D, columns=s1L, index=s2L)
    
    # Returns alingned sequences in a nice matrix format
    
    A1List = list(As1)
    A2List = list(Alg)
    A3List = list(As2)

    AlgM = np.matrix([A1List,A2List,A3List])
    
    width = 75
    i1 = [As1[i:i+width] for i in range(0, len(As1), width)]
    i2 = [Alg[i:i+width] for i in range(0, len(Alg), width)]
    i3 = [As2[i:i+width] for i in range(0, len(As2), width)]
    i4 = [Alsp[i:i+width] for i in range(0, len(Alsp), width)]
    out = [None]*(len(i1)+len(i2)+len(i3)+len(i4))
    out [::4] = i1
    out [1::4] = i2
    out [2::4] = i3
    out [3::4] = i4
    
    
    # output = ('\n'.join(map(str, out)))
            
    return M,D,As1,Alg,As2,Alsp,AlignmentScore,AlgM,out


#Results================================================================================
    
def final(Alignment, s1, s2, GapScore, ScrMtrx, SeqType, titleseq1, titleseq2):    
    print('\n')
    print("Sequence 1:")
    print('\n')
    print(titleseq1)
    print('\n')
    print("Sequence 2:")
    print('\n')
    print(titleseq2)
    print('\n')
    print("Alignment Score: " + str(Alignment(s1,s2, GapScore, ScrMtrx, SeqType)[6]))
    print('\n')
    print('\n')
    print('\n'.join(map(str, Alignment(s1,s2, GapScore, ScrMtrx, SeqType)[8])))
    
# Command line call
 
print("\n")
Atype = input("Lastly, did you want the (1) Global or (2) Local alingnment? Enter 1 or 2: ")

if Atype == '1' and SeqType == NucSeq:
    Alignment = NeedlemanWunsch
    print('\n')
    print('Global alignment coming right up. This might take a few seconds. Sit tight.')
    print('\n')
    print('\n')
    print('Match Score: ' + str(MatchScore))
    print('Mismatch SCore: ' + str(MismatchScore))
    print('Gap Score: ' + str(GapScore))
    final(Alignment, s1, s2, GapScore, ScrMtrx, SeqType, titleseq1, titleseq2)
elif Atype == '1' and SeqType == AASeq:
    Alignment = NeedlemanWunsch
    print('\n')
    print('Global alignment coming right up. This might take a few seconds. Sit tight.')
    print('\n')
    print('\n')
    print('Score Matrix: ' + str(MatrixName))
    print('Gap Score: ' + str(GapScore))
    final(Alignment, s1, s2, GapScore, ScrMtrx, SeqType, titleseq1, titleseq2)
elif Atype == '2' and SeqType == NucSeq:
    Alignment = SmithWaterman    
    print('\n')
    print('Local alignment coming right up. This might take a few seconds. Sit tight.')
    print('\n')
    print('\n')
    print('Match Score: ' + str(MatchScore))
    print('Mismatch SCore: ' + str(MismatchScore))
    print('Gap Score: ' + str(GapScore))
    final(Alignment, s1, s2, GapScore, ScrMtrx, SeqType, titleseq1, titleseq2)
elif Atype == '2' and SeqType == AASeq:
    Alignment = SmithWaterman    
    print('\n')
    print('Local alignment coming right up. This might take a few seconds. Sit tight.')
    print('\n')
    print('\n')
    print('Score Matrix: ' + str(MatrixName))
    print('Gap Score: ' + str(GapScore))
    final(Alignment, s1, s2, GapScore, ScrMtrx, SeqType, titleseq1, titleseq2)
elif SeqType == NucSeq:
    print('\n')
    print("You're pretty bad at following intructions aren't you?")
    print("Anyway, I'm going to assume you meant to say (1) Global.")
    print('This might take a few seconds. Sit tight.')
    Alignment = NeedlemanWunsch
    print('\n')
    print('\n')
    print('Match Score: ' + str(MatchScore))
    print('Mismatch SCore: ' + str(MismatchScore))
    print('Gap Score: ' + str(GapScore))
    final(Alignment, s1, s2, GapScore, ScrMtrx, SeqType, titleseq1, titleseq2)
else:
    print('\n')
    print("You're pretty bad at following intructions aren't you?")
    print("Anyway, I'm going to assume you meant to say (1) Global.")
    print('This might take a few seconds. Sit tight.')
    Alignment = NeedlemanWunsch
    print('\n')
    print('\n')
    print('Score Matrix: ' + str(MatrixName))
    print('Gap Score: ' + str(GapScore))
    final(Alignment, s1, s2, GapScore, ScrMtrx, SeqType, titleseq1, titleseq2)