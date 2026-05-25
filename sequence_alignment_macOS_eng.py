#Import the necessary libraries and modules
import sys
from Bio import SeqIO, pairwise2, Align
from Bio.Align import substitution_matrices
from Bio.pairwise2 import format_alignment

#Import your sequences from a .fasta file

seq1 = SeqIO.read('Homologous-human.fasta', 'fasta').seq
seq2 = SeqIO.read('Homologous-mouse.fasta', 'fasta').seq
shortest_sequence = min([seq1, seq2], key = len)
longest_sequence = max([seq1, seq2], key = len)

#Choose the substition matrix, uncomment the following lines by removing the '#' to see all available matrices
#print(substitution_matrices.load())
#sys.exit() 
matrix = substitution_matrices.load("PAM250")
#print(matrix)
#exit()
gap_opening_penalty = -9
gap_extension_penalty = -8.5

#Perform the pairwise alignment
alignments = pairwise2.align.localds(seq1, seq2, matrix, gap_opening_penalty, gap_extension_penalty)
my_alignment = format_alignment(*alignments[0], full_sequences=True)

#Determine the percentage of identity and similarity
lines= my_alignment.splitlines()
identities = lines[1].count('|')
similarities = lines[1].count('.')
gaps_seq1 = lines[0].count('-')
gaps_seq2 = lines[2].count('-')
id_percentage_1 = (100*identities/len(seq1))
id_percentage_2 = (100*identities/len(seq2))

#Export the results to a file
with open ('alignment.txt', 'w') as f:

        for i in range(-(-len(longest_sequence)// 90)):
                for line in lines[0:3]:
                        if 90*(i + 1) <= len(line):
                                f.write(line[90*i:90*(i+1)])
                                f.write('\n')
                        else:
                                f.write(line[90*i:-1])
                                f.write('\n')
                f.write('\n')
        f.write(lines[3].strip() + '\n')
        f.write('Sequence 1 ID%: ' + "{:.2f}".format(100*identities/len(seq1)) + '\n')
        f.write('Sequence 2 ID%: ' + "{:.2f}".format(100*lines[1].count('|')/len(seq2)) + '\n')
        f.write('Sequence 1 similarity%: ' + "{:.2f}".format(100*(lines[1].count('.')+ lines[1].count('|'))/len(seq1)) + '\n')
        f.write('Sequence 2 similarity%: ' + "{:.2f}".format(100*(lines[1].count('.')+ lines[1].count('|'))/len(seq2)) + '\n')
        f.write('Sequence 1 number of gaps: ' + str(lines[0].count('-')) + '\n')
        f.write('Sequence 2 number of gaps: ' + str(lines[2].count('-')) + '\n')


#Print a summary of the results on screen. If GREEN, the proteins are likely homologous, if RED they probably are not
if id_percentage_1 > 35:
        print('The ID% of sequence 1 is\033[1;32m',  "{:.2f}".format(id_percentage_1), '\033[00m%.')
elif id_percentage_1 >= 20 and id_percentage_1 <= 35:
        print('The ID% of sequence 1 is\033[1;33m',  "{:.2f}".format(id_percentage_1), '\033[00m%.')
elif id_percentage_1 < 20:
        print('The ID% of sequence 1 is\033[0;31m',  "{:.2f}".format(id_percentage_1), '\033[00m%.')

if id_percentage_2 > 35:
        print('The ID% of sequence 2 is\033[1;32m',  "{:.2f}".format(id_percentage_2), '\033[00m%.')
elif id_percentage_2 >= 20 and id_percentage_2 <= 35:
        print('The ID% of sequence 2 is\033[1;33m',  "{:.2f}".format(id_percentage_2), '\033[00m%.')
elif id_percentage_2 < 20:
        print('The ID% of sequence 2 is\033[0;31m',  "{:.2f}".format(id_percentage_2), '\033[00m%.')

