# Smith-Waterman=====================================================================================

# Clear console and variables if you'd like

%reset
clear()

# Import libraries

import numpy as np
import pandas as pd

# The pseudo code can be found in Wikipedia

# Dummy sequences

A = "GAATTCAGTTA"
B = "GGATCGA"

# Default Scoring System
# This can be adjusted here, or when running the 'Matrix' function (see below)

Match = 5
Mismatch = -3
Gap = -4

# The default scores for Match, Mismatch, and Gap can be adjusted when running the function

def SmithWaterman(s1,s2, MatchScore = Match, MismatchScore = Mismatch, GapScore = Gap):
    
    # Set up values for Matrices
    
    n = len(s1) # n is associated with x, i
    m = len(s2) # m is associated with y, j
    S = {True: MatchScore, False: MismatchScore} # Define S for 'Match/Mismatch' equation
    d = GapScore
    s1L = list(s1.strip()) # convert s1 string to list
    s2L = list(s2.strip()) # convert s2 string to list
    
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
    for j in range(1,len(s1)+1): # Note that here, j corresponds to s1
        for i in range(1,len(s2)+1): # Note that here, i corresponds to s2
            Match = M[i-1][j-1] + S[s1[j-1]==s2[i-1]]
            Delete = M[i-1][j] + d
            Insert = M[i][j-1] + d
            M[i,j] = max(Match, Delete, Insert,0)
            score = max(Match, Delete, Insert,0)
            if score > max_score:
                max_score = score
                max_i = [i]
                max_j = [j]
            if M[i,j] == Match:
                D[i,j] = 1 # denotes diagonal
            elif M[i,j] == Delete:
                D[i,j] = 3 # denotes vertical
            elif M[i,j] == 0:
                D[i,j] == 0
            else:
                D[i,j] = 2 # denotes horizontal
                
    # Return alignment
    
    As1 = '' # This is the stand in for alignment of s1
    As2 = '' # This is the stand in for alignment of s2
    Alg = '' # This is the stand in for alignment labels
    
    i = max_j[0]
    j = max_i[0]
    
    AlignmentScore = M[j,i]
    
    while i>0 and j>0: # Note that j comes before i
        if D[j,i] == 1:
            As1 = s1[i-1]+As1
            As2 = s2[j-1]+As2
            if s1[i-1]==s2[j-1]:
                Alg = "|"+Alg
            else:
                Alg = " "+Alg
            i = i-1
            j = j-1
        elif D[j,i] == 2:
            As1 = s1[i-1]+As1
            As2 = "-"+As2
            Alg = " "+Alg
            i = i-1
        else:
            As1 = "-"+As1
            As2 = s2[j-1]+As2
            Alg = " "+Alg
            j= j-1
    
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
            
    return M,D,As1,Alg,As2,AlignmentScore, AlgM

Seq1 = SmithWaterman(A,B)[2]
AlignmentLabels = SmithWaterman(A,B)[3]
Seq2 = SmithWaterman(A,B)[4]
AlignedSequence = SmithWaterman(A,B)[6] # Return aligned sequences in a matrix
AlignmentScore = SmithWaterman(A,B)[5]
ScoreMatrix = SmithWaterman(A,B)[0]
DirectionMatrix = SmithWaterman(A,B)[1]

