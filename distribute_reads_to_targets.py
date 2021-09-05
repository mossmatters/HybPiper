#!/usr/bin/env python

"""
This script is part of a pipeline to extract phylogenetically-useful sequences from
Illumina data using the targeted (liquid-phase) sequence enrichment approach.

After a BLASTx search of the raw reads against the target sequences, the reads need to be
sorted according to the successful hits. This script takes the BLASTx output (tabular)
and the raw read files, and distributes the reads into FASTA files ready for assembly.

If there are multiple BLAST results (for example, one for each read direction),
concatenate them prior to sorting.
"""

import sys
import os
import errno
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import argparse
import logging


# Create logger:
logger = logging.getLogger(f'__main__.{__name__}')


def mkdir_p(path):
    """

    :param path:
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

    :param blastfilename:
    :return:
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

    :param target:
    :param ID1:
    :param Seq1:
    :param Qual1:
    :param ID2:
    :param Seq2:
    :param Qual2:
    :param single:
    :param merged:
    :return:
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
    Distributing targets from single-end sequencing

    :param target:
    :param ID1:
    :param Seq1:
    :return:
    """

    mkdir_p(target)
    outfile = open(os.path.join(target, f'{target}_unpaired.fasta'), 'a')
    outfile.write(f'>{ID1}\n{Seq1}\n')
    outfile.close()


def distribute_reads(readfiles, read_hit_dict, single=True, merged=False):
    """

    :param readfiles:
    :param read_hit_dict:
    :param single:
    :param merged:
    :return:
    """

    if merged:
        logger.info('Writing fastq files for merging with BBmerge.sh')

    iterator1 = FastqGeneralIterator(open(readfiles[0]))
    if len(readfiles) == 1:

        for ID1_long, Seq1, Qual1 in iterator1:
            ID1 = ID1_long.split()[0]
            if ID1 in read_hit_dict:
                for target in read_hit_dict[ID1]:
                    write_single_seqs(target, ID1, Seq1)
        return

    elif len(readfiles) == 2:
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
    parser.add_argument('blast_filename', help='Name of blast file from BLASTX mapping of reads')
    parser.add_argument('readfiles', nargs='+', help="List of readfiles")
    parser.add_argument("--merged", help="If provided, write fastq files for bbmerge", action="store_true",
                        default=False)
    args = parser.parse_args()

    logging.info(f'Running script distribute_reads_to_targets.py with {args}')
    readfiles = args.readfiles
    read_hit_dict = read_sorting(args.blast_filename)
    logging.info(f'Unique reads with hits: {len(read_hit_dict)}')
    distribute_reads(readfiles, read_hit_dict, single=True, merged=args.merged)


if __name__ == '__main__':
    main()
