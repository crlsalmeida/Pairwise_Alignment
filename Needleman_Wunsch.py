# Needleman-Wunsch Function==================================================================================

# Clear console and variables if you'd like

%reset
clear()

# Import libraries

import numpy as np
import pandas as pd
import os

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
        
    output = print('\n'.join(map(str, out)))
            
    return M,D,As1,Alg,As2,Alsp,AlignmentScore,AlgM,out,output