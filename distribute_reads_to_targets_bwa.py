#!/usr/bin/env python

import sys, os, errno, subprocess
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

"""
This script is part of a pipeline to extract phylogenetically-useful sequences from 
Illumina data using the targeted (liquid-phase) sequence enrichment approach.

After a BWA search of the raw reads against the target sequences, the reads neeproteins had no good matchesd to be 
sorted according to the successful hits. This script takes the BWA output (BAM format)
and the raw read files, and distributes the reads into FASTA files ready for assembly.

If there are multiple results (for example, one for each read direction),
concatenate them prior to sorting.
"""


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def read_sorting(bamfilename):
    samtools_cmd = "samtools view -F 4 {}".format(bamfilename)
    child = subprocess.Popen(samtools_cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    # print(f'CJJ child: {child}')
    bwa_results = child.stdout.readlines()
    # print(f'CJJ bwa_results: {bwa_results}')

    read_hit_dict = {}
    for line in bwa_results:
        line = line.split()
        readID = line[0]
        target = line[2].split("-")[-1]
        if readID in read_hit_dict:
            if target not in read_hit_dict[readID]:
                read_hit_dict[readID].append(target)
        else:
            read_hit_dict[readID] = [target]
    return read_hit_dict


def write_paired_seqs(target, ID1, Seq1, Qual1, ID2, Seq2, Qual2, single=True):
    mkdir_p(target)
    if single:
        outfile = open(os.path.join(target, "{}_interleaved.fastq".format(target)), 'a')
        outfile.write("@{}\n{}\n{}\n{}\n".format(ID1, Seq1, f'+', Qual1))
        outfile.write("@{}\n{}\n{}\n{}\n".format(ID2, Seq2, f'+', Qual2))
        outfile.close()

        outfile = open(os.path.join(target, "{}_interleaved.fasta".format(target)), 'a')
        outfile.write(">{}\n{}\n".format(ID1, Seq1))
        outfile.write(">{}\n{}\n".format(ID2, Seq2))
        outfile.close()
    else:
        outfile1 = open(os.path.join(target, "{}_1.fasta".format(target)), 'a')
        outfile1.write(">{}\n{}\n".format(ID1, Seq1))
        outfile2 = open(os.path.join(target, "{}_2.fasta".format(target)), 'a')
        outfile2.write(">{}\n{}\n".format(ID2, Seq2))
        outfile1.close()
        outfile2.close()


def write_single_seqs(target, ID1, Seq1):
    """Distributing targets from single-end sequencing"""
    mkdir_p(target)
    outfile = open(os.path.join(target, "{}_unpaired.fasta".format(target)), 'a')
    outfile.write(">{}\n{}\n".format(ID1, Seq1))
    outfile.close()


def distribute_reads(readfiles, read_hit_dict, single=True):
    # print(f'CJJ - read_hit_dict: {read_hit_dict}')
    iterator1 = FastqGeneralIterator(open(readfiles[0]))
    if len(readfiles) == 1:

        for ID1_long, Seq1, Qual1 in iterator1:
            ID1 = ID1_long.split()[0]
            if ID1.endswith("\1") or ID1.endswith("\2"):
                ID1 = ID1[:-2]
            if ID1 in read_hit_dict:
                for target in read_hit_dict[ID1]:
                    write_single_seqs(target, ID1, Seq1)
        return

    elif len(readfiles) == 2:
        # print(f'CJJ - two reads files: { readfiles}')
        iterator2 = FastqGeneralIterator(open(readfiles[1]))

    for ID1_long, Seq1, Qual1 in iterator1:
        ID2_long, Seq2, Qual2 = next(iterator2)


        ID1 = ID1_long.split()[0]
        if ID1.endswith("/1") or ID1.endswith("/2"):
            ID1 = ID1[:-2]

        ID2 = ID2_long.split()[0]
        if ID2.endswith("/1") or ID2.endswith("/2"):
            ID2 = ID2[:-2]

        if ID1 in read_hit_dict:
            for target in read_hit_dict[ID1]:
                write_paired_seqs(target, ID1, Seq1, Qual1, ID2, Seq2, Qual2)
        elif ID2 in read_hit_dict:
            for target in read_hit_dict[ID2]:
                write_paired_seqs(target, ID1, Seq1, Qual1, ID2, Seq2, Qual2)


def main():
    bamfilename = sys.argv[1]
    # print(f'CJJ bamfilename: {bamfilename}')
    readfiles = sys.argv[2:]
    read_hit_dict = read_sorting(bamfilename)
    # print(f'CJJ:read_hit_dict {read_hit_dict}')
    print("Unique reads with hits: {}".format(len(read_hit_dict)))
    distribute_reads(readfiles, read_hit_dict, single=True)


if __name__ == "__main__": main()
