
"""AIM: The goal of this piece of code is to determine how the performance of HybPiper changes
as the sequenced reads are more diverged from the target sequence.

Tasks: We will determine this through three different steps:

1) Use Biopython to generate a series of fastq files of paired end reads where we are randomly introducing
an ever incerasing number of variation from the target sequence.

2) Test HybPiper and see how it performs. What proportions of the reads are aligned with the target and what proportion of them fail."""

from Bio import SeqIO
from Bio import Seq
import random

R1_out=open("200bases20Read_R1.fastqc","write")
R2_out=open("200bases20Read_R2.fastqc","write")

sample='@M00223:27:000000000-AAF1Y:1:1101:'
fwd_read='1:N:0:14'
rev_read='2:N:0:14'

for record in SeqIO.parse("../test_dataset/test_targets.fasta", "fasta"):
	seq = record.seq
	totLength=len(seq)
	for read in range(0,20): # We are going to simulate 20 reads for now, each=200bases long.
		n1=random.randint(20000,25000) ## this number is unique to the read
		n2=random.randint(1000,2000) ## this number is an unique to the read
		line1=sample+":"+str(n1)+":"+str(n2) ## this will be the first line of the fastq file
		point1=random.randint(0,totLength-1)
		point2=random.randint(point1+200,totLength-1)
		if point1 < (totLength - 200):
			point2=random.randint(point1+200,totLength-1)
			R=seq[point:point+200]
			line2=(str(R)
		else:
			R=seq[point-200:point]
			line2=print(str(R))	
