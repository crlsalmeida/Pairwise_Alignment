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
