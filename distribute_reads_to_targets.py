#!/usr/bin/env python

import sys, os, errno
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import argparse

"""
This script is part of a pipeline to extract phylogenetically-useful sequences from 
Illumina data using the targeted (liquid-phase) sequence enrichment approach.

After a BLASTx search of the raw reads against the target sequences, the reads need to be 
sorted according to the successful hits. This script takes the BLASTx output (tabular)
and the raw read files, and distributes the reads into FASTA files ready for assembly.

If there are multiple BLAST results (for example, one for each read direction),
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


def read_sorting(blastfilename):
    read_hit_dict = {}
    blastfile = open(blastfilename)
    for line in blastfile:
        line = line.split()
        readID = line[0]
        target = line[1].split("-")[-1]
        if readID in read_hit_dict:
            if target not in read_hit_dict[readID]:
                read_hit_dict[readID].append(target)
        else:
            read_hit_dict[readID] = [target]
    return read_hit_dict


def write_paired_seqs(target, ID1, Seq1, Qual1, ID2, Seq2, Qual2, single=True, merged=False):
    mkdir_p(target)
    if single:
        if merged:
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
    """
    Distributing targets from single-end sequencing
    """

    mkdir_p(target)
    outfile = open(os.path.join(target, "{}_unpaired.fasta".format(target)), 'a')
    outfile.write(">{}\n{}\n".format(ID1, Seq1))
    outfile.close()


def distribute_reads(readfiles, read_hit_dict, single=True, merged=False):
    iterator1 = FastqGeneralIterator(open(readfiles[0]))
    if len(readfiles) == 1:

        for ID1_long, Seq1, Qual1 in iterator1:
            ID1 = ID1_long.split()[0]
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
        ID2 = ID2_long.split()[0]

        if ID1 in read_hit_dict:
            for target in read_hit_dict[ID1]:
                write_paired_seqs(target, ID1, Seq1, Qual1, ID2, Seq2, Qual2, merged=merged)
        elif ID2 in read_hit_dict:
            for target in read_hit_dict[ID2]:
                write_paired_seqs(target, ID1, Seq1, Qual1, ID2, Seq2, Qual2, merged=merged)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("blast_filename", help="Name of blast file from BLASTX mapping of reads")  # CJJ
    parser.add_argument('readfiles', nargs='+', help="List of readfiles")  # CJJ
    parser.add_argument("--merged", help="If provided, write fastq files for bbmerge",
                        action="store_true", default=False)  # CJJ
    args = parser.parse_args()

    # blastfilename = sys.argv[1]
    # readfiles = sys.argv[2:]
    readfiles = args.readfiles
    read_hit_dict = read_sorting(args.blast_filename)
    # print read_hit_dict
    print("Unique reads with hits: {}".format(len(read_hit_dict)))
    distribute_reads(readfiles, read_hit_dict, single=True, merged=args.merged)


if __name__ == "__main__":
    main()
