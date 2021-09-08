#!/usr/bin/env python

"""
This script is part of a pipeline to extract phylogenetically-useful sequences from
Illumina data using the targeted (liquid-phase) sequence enrichment approach.

After a BLASTx search of the raw reads against the target sequences, the reads need to be
sorted according to the successful hits. This script takes the BLASTx output (tabular)
and the raw read files, and distributes the reads into FASTA files ready for assembly.

If there are multiple BLAST results (for example, one for each read direction),
concatenate them prior to sorting. # CJJ not still true?
"""

import os
import errno
import argparse
import logging
from distribute_reads_to_targets_bwa import distribute_reads


# Create logger:
logger = logging.getLogger(f'__main__.{__name__}')


def mkdir_p(path):
    """
    Creates a directory corresponding the the given path, if it doesn't already exist.

    :param str path: path of directory to create
    :return:
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def read_sorting(blastfilename):
    """
    Returns a dictionary of read_hit_dict[readID] = [target1, target2, ...]

    :param str blastfilename: path the BLASTx tabular output file
    :return: dict read_hit_dict: dictionary of read_hit_dict[readID] = [target1, target2, ...]
    """

    read_hit_dict = {}
    blastfile = open(blastfilename)
    for line in blastfile:
        line = line.split()
        readID = line[0]
        target = line[1].split('-')[-1]
        if readID in read_hit_dict:
            if target not in read_hit_dict[readID]:
                read_hit_dict[readID].append(target)
        else:
            read_hit_dict[readID] = [target]
    return read_hit_dict


def write_paired_seqs(target, ID1, Seq1, Qual1, ID2, Seq2, Qual2, single=True, merged=False):
    """
    Writes interleaved fasta files, and also interleaved fastq files if merged=True, to the corresponding gene directory

    :param str target: gene name e.g. gene001
    :param str ID1: fasta/fastq header for R1
    :param str Seq1: fasta/fastq sequence for R1
    :param str Qual1: fastq quality scores for R1
    :param str ID2: fasta/fastq head for R2
    :param str Seq2: fasta/fastq sequence for R2
    :param str Qual2: fastq quality scores for R2
    :param bool single: # CJJ hardcoded as True - remove?
    :param bool merged: If True, write fastq seqs as well as fasta
    :return::
    """

    mkdir_p(target)
    if single:  # If True, write paired reads in interleaved format
        if merged:
            outfile = open(os.path.join(target, f'{target}_interleaved.fastq'), 'a')
            outfile.write(f'@{ID1}\n{Seq1}\n+\n{Qual1}\n')
            outfile.write(f'@{ID2}\n{Seq2}\n+\n{Qual2}\n')
            outfile.close()

        outfile = open(os.path.join(target, f'{target}_interleaved.fasta'), 'a')
        outfile.write(f'>{ID1}\n{Seq1}\n')
        outfile.write(f'>{ID2}\n{Seq2}\n')
        outfile.close()
    # else:
    #     outfile1 = open(os.path.join(target, "{}_1.fasta".format(target)), 'a')
    #     outfile1.write(">{}\n{}\n".format(ID1, Seq1))
    #     outfile2 = open(os.path.join(target, "{}_2.fasta".format(target)), 'a')
    #     outfile2.write(">{}\n{}\n".format(ID2, Seq2))
    #     outfile1.close()
    #     outfile2.close()


def write_single_seqs(target, ID1, Seq1):
    """
    Writes a fasta filesfo single end reads to the corresponding gene directory

    :param str target: gene name e.g. gene001
    :param str ID1: fasta/fastq header for R1
    :param str Seq1: fasta/fastq sequence for R1
    :return:
    """

    mkdir_p(target)
    outfile = open(os.path.join(target, f'{target}_unpaired.fasta'), 'a')
    outfile.write(f'>{ID1}\n{Seq1}\n')
    outfile.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('blast_filename', help='Name of blast file from BLASTX mapping of reads')
    parser.add_argument('readfiles', nargs='+', help="List of readfiles")
    parser.add_argument("--merged", help="If provided, write fastq files for bbmerge", action="store_true",
                        default=False)
    args = parser.parse_args()

    logging.info(f'Running script distribute_reads_to_targets.py with {args}')
    readfiles = args.readfiles
    read_hit_dict = read_sorting(args.blast_filename)
    logging.info(f'Unique reads with hits: {len(read_hit_dict)}')
    distribute_reads(readfiles, read_hit_dict, merged=args.merged)


if __name__ == '__main__':
    main()
