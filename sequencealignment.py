"""

Align a DNA sequence and a protein sequence based on the translated sequences of the DNA sequences.
Only the forward reading frames need be considered.
If desired, you may use existing packages for reading/writing sequences, scoring tables, and formatting the results.
Display the results in conventional pairs of lines format, i.e. in blocks with the aligned sequences above each other, showing both the DNA and translated protein sequences.
 Use a local alignment algorithm
o Bonus: +10 pts. Implement both global and local alignments
 Assume the sequences are in fasta format
 Use a scoring table such as Blosum or PAM read in from a file. A good source of scoring
matrix files is ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/
 Allow affine gaps within reading frames
 Allow length independent gaps between reading frames.
o Bonus + 5 pts. : allow affine gaps between reading frames

"""
import sys
from Bio import SeqIO
import numpy as np
from numpy import unravel_index


dnafasta=sys.argv[1]
proteinfasta=sys.argv[2]
flag=sys.argv[3]
penalty=float(sys.argv[4])
scoringmatrix=sys.argv[5]

for seq_record in SeqIO.parse(dnafasta, "fasta"):
  str1=str(seq_record.seq.translate())
  #print(str1)

for seq_record in SeqIO.parse(proteinfasta, "fasta"):
  str2=str(seq_record.seq)

seq1=str1.upper()
seq2=str2.upper()

class Matrix:
  def __init__(self, matrix_filename):
    self._load_matrix(matrix_filename)

  def _load_matrix(self, matrix_filename):
    with open(matrix_filename) as matrix_file:
      matrix = matrix_file.read()
    lines = matrix.strip().split('\n')

    header = lines.pop(0)
    columns = header.split()
    matrix = {}

    for row in lines:
      entries = row.split()
      row_name = entries.pop(0)
      matrix[row_name] = {}

      if len(entries) != len(columns):
        raise Exception('Improper entry number in row')
      for column_name in columns:
        matrix[row_name][column_name] = entries.pop(0)

    self._matrix = matrix

  def lookup_score(self, a, b):
        a = a.upper()
        b = b.upper()

        if a not in self._matrix or b not in self._matrix[a]:
            raise InvalidPairException('[%s, %s]' % (a, b))
        return self._matrix[a][b]

BlosumMatrix=Matrix(scoringmatrix)
BlosumMatrix._load_matrix(scoringmatrix)


print("Length of translated protein sequence is "+str(len(seq1)))
print("Length of given protein sequence is "+str(len(seq2)))

""""

seq1="PQRLGSSAAAAAAAAAAAAAQRDSTYQQQQ"
seq2="AAAAAAAAAAQR*S*S*RL**LGGTD"
"""

#print(seq1)
#print(seq2)
rows, cols = (len(seq2)+2,len(seq1)+2)

matrix =[[0]*cols for _ in range(rows)]

if(int(flag)==1):
    for y in range(2, len(seq2) + 2):
        matrix[y][0] = (seq2[y - 2])
        matrix[y][1] = penalty * (y - 1)

    for y in range(2, len(seq1) + 2):
        matrix[0][y] = (seq1[y - 2])
        matrix[1][y] = penalty * (y - 1)

    matrix[0][1] = '_'
    matrix[1][0] = '_'

    rows, cols = (len(seq2) + 2, len(seq1) + 2)
    traceback = [[0] * cols for _ in range(rows)]

    for y in range(2, len(seq2) + 2):
        traceback[y][0] = (seq2[y - 2])

    for y in range(2, len(seq1) + 2):
        traceback[0][y] = (seq1[y - 2])

    traceback[0][1] = '_'
    traceback[1][0] = '_'

    print("Global Alignment")

    for i in range(2, len(seq2) + 2):
        for j in range(2, len(seq1) + 2):

            value = float(BlosumMatrix.lookup_score(matrix[i][0], matrix[0][j]))
            notgap = matrix[i - 1][j - 1] + value
            gap1 = matrix[i - 1][j] + penalty
            gap2 = matrix[i][j - 1] + penalty
            final = max(notgap, gap1, gap2)

            if final==notgap:
                if (matrix[i][0] == matrix[0][j]):
                    traceback[i][j]=("M->D")
                else:
                    traceback[i][j]=("MM->D")
            elif final==gap1:
                traceback[i][j]=("G->T")
            elif final == gap2:
                traceback[i][j]=("G->L")

            matrix[i][j] = final


    print(np.matrix(matrix))
    alignedseq1=[]
    alignedseq2=[]

    print("Global Alignment score: "+str(matrix[len(seq2)+1][len(seq1)+1]))

    while(i>=1 and j>=1):

            if((traceback[i][j]=="MM->D") | (traceback[i][j]=="M->D") ):

               alignedseq2.append(traceback[i][0])
               alignedseq1.append(traceback[0][j])
               j=j-1
               i=i-1

            elif((traceback[i][j]==0)):
                if(i==1 and j==1):
                    break

                elif(i==1):
                    alignedseq2.append(traceback[i][0])
                    alignedseq1.append(traceback[0][j])
                    j=j-1
                elif(j==1):
                    alignedseq2.append(traceback[i][0])
                    alignedseq1.append(traceback[0][j])
                    i=i-1

            elif((traceback[i][j]=="G->T")):
               alignedseq1.append('_')
               alignedseq2.append(traceback[i][0])
               i=i-1
            elif((traceback[i][j]=="G->L")):
               alignedseq2.append('_')
               alignedseq1.append(traceback[0][j])
               j=j-1


    print(alignedseq1[::-1])
    print(alignedseq2[::-1])


elif(int(flag)==2):
    print("Local Alignment")
    for y in range(2, len(seq2) + 2):
        matrix[y][0] = (seq2[y - 2])
        matrix[y][1] = 0

    for y in range(2, len(seq1) + 2):
        matrix[0][y] = (seq1[y - 2])
        matrix[1][y] = 0

    matrix[0][1] = '_'
    matrix[1][0] = '_'

    rows, cols = (len(seq2) + 2, len(seq1) + 2)
    traceback = [[0] * cols for _ in range(rows)]

    for y in range(2, len(seq2) + 2):
        traceback[y][0] = (seq2[y - 2])

    for y in range(2, len(seq1) + 2):
        traceback[0][y] = (seq1[y - 2])

    traceback[0][1] = '_'
    traceback[1][0] = '_'

    for i in range(2, len(seq2) + 2):
        for j in range(2, len(seq1) + 2):

            value=float(BlosumMatrix.lookup_score(matrix[i][0],matrix[0][j]))

            notgap = matrix[i - 1][j - 1] + value
            gap1 = matrix[i - 1][j] + penalty
            gap2 = matrix[i][j - 1] + penalty
            final = max(notgap, gap1, gap2,0)

            if final==notgap:
                if (matrix[i][0] == matrix[0][j]):
                    traceback[i][j]=("M->D")
                else:
                    traceback[i][j]=("MM->D")
            elif final==gap1:
                traceback[i][j]=("G->T")
            elif final == gap2:
                traceback[i][j]=("G->L")

            matrix[i][j] = final


    print(np.matrix(matrix))
    alignedseq1=[]
    alignedseq2=[]

    matrix = np.array(matrix)

    dummymatrix = np.delete(matrix, 0, 0)
    dummymatrix = np.delete(dummymatrix, 0,1 )

    max_tuple=unravel_index(dummymatrix.argmax(), dummymatrix.shape)

    i=max_tuple[0]+1
    j=max_tuple[1]+1
    print("Local Alignment score: " + str(dummymatrix[max_tuple[0]][max_tuple[1]]))
    while( j>=1 and i>=1):

            if((traceback[i][j]=="MM->D") | (traceback[i][j]=="M->D") ):

               alignedseq2.append(traceback[i][0])
               alignedseq1.append(traceback[0][j])
               j=j-1
               i=i-1

            elif((traceback[i][j]==0) ):
                break

            elif((traceback[i][j]=="G->T")):
               alignedseq1.append('_')
               alignedseq2.append(traceback[i][0])
               i=i-1
            elif((traceback[i][j]=="G->L")):
               alignedseq2.append('_')
               alignedseq1.append(traceback[0][j])
               j=j-1


    print(alignedseq1[::-1])
    print(alignedseq2[::-1])








