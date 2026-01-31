"""AIM: The goal of this script is to randomly introduce variation into a series of fastq files.
"""

from Bio import SeqIO
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq

inFile = SeqIO.parse("/Users/laurabotigue/HybPiper/test_dataset/EG30_R1_test.fastq", "fastq")
outFile= open("/Users/laurabotigue/HybPiper/test_dataset/EG30_R1_iter2.fastq", "w")

bases=['A','C','G','T']

for read in inFile:
	len_seq=len(read.seq)
	pos = random.randint(0,len_seq-1)
	bases.remove(str(read.seq[pos])) #check if using dictionary instead of list is more elegant.
	mut = random.sample(bases, 1) 
	output_seq = read.seq.tomutable()
	output_seq[pos] = mut[0]
	bases=['A','C','G','T'] ## This is to restore the poss number of mutations

	new_seq = SeqRecord(output_seq, id=read.id, description=read.description)
	new_seq.letter_annotations["phred_quality"]=read.letter_annotations["phred_quality"] # A rather convoluted way to pair the new sequence with the old fastq information	
	
	SeqIO.write(new_seq, outFile, "fastq")
	
outFile.close()


