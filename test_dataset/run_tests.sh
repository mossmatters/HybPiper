#!/bin/bash


#Detect if test_reads are the real thing or git-lfs placeholder
minsize=5000
readsize=$(wc -c < test_reads.fastq.tar.gz)
if [ $minsize -ge $readsize ]; then
    rm test_reads.fastq.tar.gz
    wget https://github.com/mossmatters/HybPiper/raw/develop/test_dataset/test_reads.fastq.tar.gz || curl -O https://github.com/mossmatters/HybPiper/raw/develop/test_dataset/test_reads.fastq.tar.gz
fi

#Unpack the test dataset
tar -zxf test_reads.fastq.tar.gz

#Remove any previous runs
parallel rm -r {} :::: namelist.txt


#Run main HybPiper script with all available CPUs
while read i
do
../reads_first.py -r $i*.fastq -b test_targets.fasta --prefix $i --bwa 
done < namelist.txt

#Get the seq_lengths.txt file
python ../get_seq_lengths.py test_targets.fasta namelist.txt dna > test_seq_lengths.txt

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
