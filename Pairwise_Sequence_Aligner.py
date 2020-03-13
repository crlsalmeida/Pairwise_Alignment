
"""
Created on Fri Jan 31 11:16:25 2020
Final edits on March 10 2020

@author: crlsalmeida
"""

# Analysis Preparation========================================================================================

# Clear console and variables if you'd like

%reset
clear()

# Import libraries

import numpy as np
import pandas as pd
import os

# Define Working Directory
# Example: "C:\\Users\\crlsalmeida\\Documents\\Programing Projects"
# Example: "C:\\Users\\Carlos\\Documents\\School and Work\\UofL\\Courses\\UofL Spring 2020\\CECS660\\Programing Projects"
workingdirectory = "C:\\Users\\Carlos\\Documents\\School and Work\\UofL\\Courses\\UofL Spring 2020\\CECS660\\Programing Projects"
os.chdir(workingdirectory)

# Check that Working Directory matches "workingdirectory"
# pwd

# Enter the names of the 2 FASTA files for comparison
# Example: DNA_seq1.FASTA, DNA_seq2.FASTA
fasta1 = "AA_seq1.FASTA"
fasta2 = "AA_seq2.FASTA"

# Run the following section to retrive the sequence from FASTA files as strings as variables A and B

def read_fasta(fasta1, fasta2, cwd = workingdirectory):
    
    '''
    This function opens 2 FASTA files
    It will return the first line of both FASTA files as 'title1' and 'title2'
    It will return both sequences as strings
    
    '''
    os.chdir(workingdirectory)
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
       
A = read_fasta(fasta1, fasta2)[0]
B = read_fasta(fasta1, fasta2)[1]
titleseq1 = read_fasta(fasta1, fasta2)[2]
titleseq2 = read_fasta(fasta1, fasta2)[3]

# At this point you may delete unnecessary items from your environment
del(fasta1, fasta2)

# Score Matrices and Look-up Functions==================================================================================

# Define Gap Score

GapScore = -4

# Default Nucleotide Scoring values and table

MatchScore = 5
MismatchScore = -3

NucScrTbl = np.array([[MatchScore, MismatchScore, MismatchScore, MismatchScore],
                      [MismatchScore, MatchScore, MismatchScore, MismatchScore],
                      [MismatchScore, MismatchScore, MatchScore, MismatchScore],
                      [MismatchScore, MismatchScore, MismatchScore, MatchScore]], dtype = float)

# Amino acid score matrices
    
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


# The following functions will be used in the Needleman-Wunsch and Smith-Waterman Functions
    
# Below is the function used to aligning nucleotide sequences
    
def NucSeq (ch1, ch2, ScrMtrx = NucScrTbl):
    
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

# Below is the function used to alignment of amino acid sequences

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

# Needleman-Wunsch Function==================================================================================

# ScrMtrx corresponds to the Score Matrix for the alignment
# By default, this function is set to align nucleotide sequences
# To align sequencues, specify the ScrMtrx (PAM250, BLOSUM62, NucScrTbl) and SeqType (AASeq or NucSeq)

def NeedlemanWunsch(s1,s2, GapScore = GapScore, ScrMtrx = NucScrTbl, SeqType = NucSeq):
    
    '''
    
    This function performs a Needleman-Wunsch global sequence alignment without affine gap penalty
    By default, this is set up to align nucleotide sequences, SeqType
        Define the SeqType as either 'NucSeq' for nucleotide sequences
        Or 'AASeq' for amino acid sequences
    The default scoring matrix, ScrMtrx, is a nucleotide score matrix
        Define ScrMtrx as either 'NucScrTbl' for nucleotide sequences
        Or 'PAM250', 'BLOSUM62', etc. for AA sequences

    '''
    
    # Set up values for Matrices
    
    n = len(s1) # n is associated with x, i
    m = len(s2) # m is associated with y, j
    s1L = list(s1.strip()) # convert s1 string to list
    s2L = list(s2.strip()) # convert s2 string to list
    d = GapScore
    
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
    
    # Returns alingned sequences
    
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
            
    return M,D,As1,Alg,As2,Alsp,AlignmentScore,AlgM,out

# Smith-Waterman Function==================================================================================

# ScrMtrx corresponds to the Score Matrix for the alignment
# By default, this function is set to align nucleotide sequences
# To align sequencues, specify the ScrMtrx (PAM250, BLOSUM62, NucScrTbl) and SeqType (AASeq or NucSeq) 

def SmithWaterman(s1,s2, GapScore = GapScore, ScrMtrx = NucScrTbl, SeqType = NucSeq):
    
    '''
    
    This function performs a Needleman-Wunsch global sequence alignment without affine gap penalty
    By default, this is set up to align nucleotide sequences, SeqType
        Define the SeqType as either 'NucSeq' for nucleotide sequences
        Or 'AASeq' for amino acid sequences
    The default scoring matrix, ScrMtrx, is a nucleotide score matrix
        Define ScrMtrx as either 'NucScrTbl' for nucleotide sequences
        Or 'PAM250', 'BLOSUM62', etc. for AA sequences

    '''
    
    # Set up values for Matrices
    
    n = len(s1) # n is associated with x, i
    m = len(s2) # m is associated with y, j
    s1L = list(s1.strip()) # convert s1 string to list
    s2L = list(s2.strip()) # convert s2 string to list
    d = GapScore

        
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
    
    # Returns alingned sequences
    
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
            
    return M,D,As1,Alg,As2,Alsp,AlignmentScore,AlgM,out

# Results=======================================================================================================

#Return Nucleuotide sequences

# (s1,s2, GapScore = GapScore, ScrMtrx = NucScrTbl, SeqType = NucSeq)

Local_Seq1Aligned = SmithWaterman(A,B)[2]
Local_Seq2Aligned = SmithWaterman(A,B)[4]
Local_Seq12Symbols = SmithWaterman(A,B)[3]
Local_SeqSpace = SmithWaterman(A,B)[5]
Local_ScoreMatrix = SmithWaterman(A,B)[0]
Local_TracebackMatrix = SmithWaterman(A,B)[1]
Local_AlignmentScore = SmithWaterman(A,B)[6]
Local_SeqAlignedTable = SmithWaterman(A,B)[7]
Local_SeqAlignedArray = SmithWaterman(A,B)[8]
print('\n'.join(map(str, SmithWaterman(A,B)[8])))

Glocal_Seq1Aligned = NeedlemanWunsch(A,B)[2]
Global_Seq2Aligned = NeedlemanWunsch(A,B)[4]
Global_Seq12Symbols = NeedlemanWunsch(A,B)[3]
Global_SeqSpace = NeedlemanWunsch(A,B)[5]
Global_ScoreMatrix = NeedlemanWunsch(A,B)[0]
Global_TracebackMatrix = NeedlemanWunsch(A,B)[1]
Global_AlignmentScore = NeedlemanWunsch(A,B)[6]
Global_SeqAlignedTable = NeedlemanWunsch(A,B)[7]
Global_SeqAlignedArray = NeedlemanWunsch(A,B)[8]
print('\n'.join(map(str, SmithWaterman(A,B)[8])))

#Return Amino Acids sequences

# (s1,s2, GapScore = GapScore, ScrMtrx = BLOSUM80(or PAM250 or BLOSUM62), SeqType = AASeq)

Local_Seq1Aligned = SmithWaterman(A,B, GapScore = GapScore, ScrMtrx = BLOSUM80, SeqType = AASeq)[2]
Local_Seq2Aligned = SmithWaterman(A,B, GapScore = GapScore, ScrMtrx = BLOSUM80, SeqType = AASeq)[4]
Local_Seq12Symbols = SmithWaterman(A,B, GapScore = GapScore, ScrMtrx = BLOSUM80, SeqType = AASeq)[3]
Local_SeqSpace = SmithWaterman(A,B, GapScore = GapScore, ScrMtrx = BLOSUM80, SeqType = AASeq)[5]
Local_ScoreMatrix = SmithWaterman(A,B, GapScore = GapScore, ScrMtrx = BLOSUM80, SeqType = AASeq)[0]
Local_TracebackMatrix = SmithWaterman(A,B, GapScore = GapScore, ScrMtrx = BLOSUM80, SeqType = AASeq)[1]
Local_AlignmentScore = SmithWaterman(A,B, GapScore = GapScore, ScrMtrx = BLOSUM80, SeqType = AASeq)[6]
Local_SeqAlignedArray = SmithWaterman(A,B, GapScore = GapScore, ScrMtrx = BLOSUM80, SeqType = AASeq)[7]

Glocal_Seq1Aligned = NeedlemanWunsch(A,B, GapScore = GapScore, ScrMtrx = BLOSUM80, SeqType = AASeq)[2]
Global_Seq2Aligned = NeedlemanWunsch(A,B, GapScore = GapScore, ScrMtrx = BLOSUM80, SeqType = AASeq)[4]
Global_Seq12Symbols = NeedlemanWunsch(A,B, GapScore = GapScore, ScrMtrx = BLOSUM80, SeqType = AASeq)[3]
Global_SeqSpace = NeedlemanWunsch(A,B, GapScore = GapScore, ScrMtrx = BLOSUM80, SeqType = AASeq)[5]
Global_ScoreMatrix = NeedlemanWunsch(A,B, GapScore = GapScore, ScrMtrx = BLOSUM80, SeqType = AASeq)[0]
Global_TracebackMatrix = NeedlemanWunsch(A,B, GapScore = GapScore, ScrMtrx = BLOSUM80, SeqType = AASeq)[1]
Global_AlignmentScore = NeedlemanWunsch(A,B, GapScore = GapScore, ScrMtrx = BLOSUM80, SeqType = AASeq)[6]
Global_SeqAlignedArray = NeedlemanWunsch(A,B, GapScore = GapScore, ScrMtrx = BLOSUM80, SeqType = AASeq)[7]

# Save to text file

# Name the text file
# ex: "MyAlignedSequences.txt"
title = 'MyAlignedSequences.txt'

# Alignment refers to the type of alignment, NeedlemanWunsh or SmithWaterman

def save_output (title, Alignment, s1, s2, GapScore, ScrMtrx, SeqType, titleseq1, titleseq2):
    output = open(title,'w+')
    output.write('\n')
    output.write("Sequence 1:")
    output.write('\n')
    output.write(titleseq1)
    output.write('\n')
    output.write("Sequence 2:")
    output.write('\n')
    output.write(titleseq2)
    output.write('\n')
    output.write('\n')
    output.write("Alignment Score:")
    output.write('\n')
    output.write(str(Alignment(s1,s2, GapScore, ScrMtrx, SeqType)[6]))
    output.write('\n')
    output.write('\n')
    output.write('\n'.join(map(str, Alignment(s1,s2, GapScore, ScrMtrx, SeqType)[8])))
    output.close()
    print("Task completed.")
 
    

# save_output (title, Alignment, s1, s2, GapScore, ScrMtrx, SeqType, titleseq1, titleseq2)

