1. Alignments DNA+Protein alignment (40 pts-55 pts)

Align a DNA sequence and a protein sequence based on the translated sequences of the DNA sequences. Only the forward reading frames need be considered.
If desired, you may use existing packages for reading/writing sequences, scoring tables, and formatting the results.
Display the results in conventional pairs of lines format, i.e. in blocks with the aligned sequences above each other, showing both the DNA and translated protein sequences.



Approach-

For this task, I used Python3 . I used NWMod to code this assignment. 

How to run the program: python3  sequencealignment.py filed.fasta filep.fasta 1 -5 blosum62.txt

In above example, the name of my program is sequencealignment.py. file1.fasta houses the DNA sequence that gets translated. file2.fasta houses the protein sequence. 1 is the flag for Global Alignment( 2 is used for Local alignment). -5 is the gap penalty. blosum62.txt is the text file from which blosum scores are read.



I have used a length independent gap penalty here.

Apart from that, I am using Python libraries like BioPython and NumPy for reading fasta files and matrix manipulation.

My code prints out the length of the translated and given protein sequence. Next, I try printing the DP table. I also report the global alignment score( or local alignment score). Finally, I align both the sequences together.



NOTE: A BioPython warning may show up during translation. This warning can be ignored.


