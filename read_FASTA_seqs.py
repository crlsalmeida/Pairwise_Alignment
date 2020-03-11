# Return Sequence from FASTA File=====================================================================================

# Clear console and variables if you'd like

%reset
clear()

# Import libraries

import os


# Define Working Directory
# Example: "C:\\Users\\Carlos\\Documents\\Folder1\\Folder2\\Folder3"
workingdirectory = "¯\_(ツ)_/¯"
os.chdir(workingdirectory)

# Check that Working Directory matches "workingdirectory"
pwd

# Enter the names of the 2 FASTA files for comparison
# Example: "filename1.fasta", "filename2.fasta"
fasta1 = "XOXOXO.fasta"
fasta2 = "YOYOYO.fasta"

# Run the following section to retrive the sequence from fasta files as strings variables A and B

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
