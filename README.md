# Pairwise_Alignment
Needleman-Wunsh and Smith-Waterman scripts in Python

Author Notes:

The scripts in this folder were part of a graduate level introductory bioinformatics project.

Feel free to use my codes, build upon them, etc.

For those who can read, understand, and execute my code, you will get the correct alignment and aligment scores at the very least.

Note, there are different ways to traceback the alignment. Depending on how one writes the code and prioritizes the traceback,
the alignment and gap locations may vary. That is, my traceback (and alignment) may vary from somebody else's at a few locations.

Files:

Cmd_promp_Pairwise_Sequence_Aligner.py
  - This file allows the user to align 2 sequences in FASTA format using a command prompt window
 
Pairwise_Sequence_Aligner.py
  - This file is the complete code for the Pairwise aligner, and contains both Needleman Wunsch and Smith Waterman alignment codes  

Needleman_Wunsch.py
  - This file only contains the code for the Global Alignment

Smith_Waterman.py
  - This file only contains the code the Local Alignment
  
 Results.py
   - This file contains the code from the 'Results' section of Pairwise_Sequence_Aligner.py
 
 Score_Matrices_&_Look-up_Functions.py
   - This file contains the code from the 'Score Matrices and Look-up Functions' section of Pairwise_Sequence_Aligner.py
 
 read_FASTA_seqs.py
   - This file contains the code from the 'Analysis Preparation' secion of Pairwise_Sequence_Aligner.py
 
 DNA_seq1.FASTA, DNA_seq2.FASTA, AA_seq1.FASTA, AA_seq2.FASTA
   - These are FASTA files downloaded from NCBI placed here to allow you test my code
   - For DNA sequences, Match=5 Mismatch=-3 Gap=4, Local Score is 1413 and Global Score is 500
   - For the AA sequences, BLOSUM80 Gap=-4, Local Score is 233 and Global Score is 203
   
