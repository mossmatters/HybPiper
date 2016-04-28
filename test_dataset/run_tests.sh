#!/bin/bash

#Unpack the test dataset
tar -zxf test_reads.fastq.tar.gz

#Run main HybPiper script with all available CPUs
while read i
do
../reads_first.py -r "$i*.fastq" -b test_targets.fasta --prefix $i --bwa 
done < namelist.txt

#Get the seq_lengths.txt file
python ../get_seq_lengths.py namelist.txt test_targets.fasta dna > test_seq_lengths.txt

#Test for paralogs
while read i
do
python ../paralog_investigator.py $i
done < namelist.txt

#Run intronerate
while read i
do
python ../intronerate.py --prefix $i
done < namelist.txt
