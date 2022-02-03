#!/usr/bin/env python

"""
This script is part of a pipeline to extract phylogenetically-useful sequences from 
Illumina data using the targeted (liquid-phase) sequence enrichment approach.

After a BWA search of the raw reads against the target sequences, the reads need to be 
sorted according to the successful hits. This script takes the BWA output (BAM format)
and the raw read files, and distributes the reads into FASTA files ready for assembly.

If there are multiple results (for example, one for each read direction),
concatenate them prior to sorting.
"""

import sys
import os
import errno
import subprocess
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import argparse
import logging
import gzip
import progressbar


# Create logger:
logger = logging.getLogger(f'hybpiper.hybpiper_assemble.{__name__}')


def mkdir_p(path):
    """
    Creates a directory corresponding to the given path, if it doesn't already exist.

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


def read_sorting(bamfilename):
    """
    Returns a dictionary of read_hit_dict[readID] = [target1, target2, ...]

    :param str bamfilename: path the *.bam file output by BWA
    :return: dict read_hit_dict: dictionary of read_hit_dict[readID] = [target1, target2, ...]
    """

    samtools_cmd = f'samtools view -F 4 {bamfilename}'
    child = subprocess.Popen(samtools_cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    bwa_results = child.stdout.readlines()

    read_hit_dict = {}
    for line in bwa_results:
        line = line.split()
        readID = line[0]  # e.g. A00119:385:H3GYVDSX2:2:1101:1416:1000
        target = line[2].split('-')[-1]  # e.g. 6128
        if readID in read_hit_dict:
            if target not in read_hit_dict[readID]:
                read_hit_dict[readID].append(target)
                # i.e. could have key(A00119:385:H3GYVDSX2:2:1101:1416:1000):value(6128, 5968)
        else:
            read_hit_dict[readID] = [target]
    return read_hit_dict


def write_paired_seqs(target, ID1, Seq1, Qual1, ID2, Seq2, Qual2, single=True, merged=False):
    """
    Writes interleaved fasta files, and also interleaved fastq files if merged=True

    :param str target: gene name e.g. gene001
    :param str ID1: fasta/fastq header for R1
    :param str Seq1: fasta/fastq sequence for R1
    :param str Qual1: fastq quality scores for R1
    :param str ID2: fasta/fastq head for R2
    :param str Seq2: fasta/fastq sequence for R2
    :param str Qual2: fastq quality scores for R2
    :param bool single: # CJJ hardcoded as True - remove?
    :param bool merged: If True, write fastq seqs as well as fasta
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


def write_single_seqs(target, ID1, Seq1):
    """
    Writes a fasta files of single-end reads to the corresponding gene directory

    :param str target: gene name e.g. gene001
    :param str ID1: fasta/fastq header for R1
    :param str Seq1: fasta/fastq sequence for R1
    :return::
    """

    mkdir_p(target)
    outfile = open(os.path.join(target, f'{target}_unpaired.fasta'), 'a')
    outfile.write(f'>{ID1}\n{Seq1}\n')
    outfile.close()
    
    
def distribute_reads(readfiles, read_hit_dict, merged=False):
    """

    :param list readfiles: a list of one or more readfiles
    :param dict read_hit_dict: dictionary of read_hit_dict[readID] = [target1, target2, ...]
    :param bool merged: boolean passed to function write_paired_seqs()
    :return:
    """

    if merged:
        logger.info(f'{"[NOTE]:":10} Writing fastq files for merging with BBmerge.sh')

    num_reads_in_readfile = 0

    # Check if read file is gzipped:
    filename, file_extension = os.path.splitext(readfiles[0])
    if file_extension == '.gz':
        logger.debug(f'Distributing reads from gzipped file {os.path.basename(readfiles[0])}')
        iterator1 = FastqGeneralIterator(gzip.open(readfiles[0], 'rt'))
        for read in iterator1:
            num_reads_in_readfile += 1
        iterator1 = FastqGeneralIterator(gzip.open(readfiles[0], 'rt'))

    else:
        iterator1 = FastqGeneralIterator(open(readfiles[0]))
        for read in iterator1:
            num_reads_in_readfile += 1
        iterator1 = FastqGeneralIterator(open(readfiles[0]))

    if len(readfiles) == 1:
        logger.info(f'{"[NOTE]:":10} Distributing unpaired reads to gene directories')

        for ID1_long, Seq1, Qual1 in progressbar.progressbar(iterator1, max_value=num_reads_in_readfile,
                                                             min_poll_interval=30):
            ID1 = ID1_long.split()[0]
            if ID1.endswith('\1') or ID1.endswith('\2'):
                ID1 = ID1[:-2]
            if ID1 in read_hit_dict:
                for target in read_hit_dict[ID1]:
                    write_single_seqs(target, ID1, Seq1)
        return

    elif len(readfiles) == 2:
        logger.info(f'{"[NOTE]:":10} Distributing paired reads to gene directories')

        # Check if read file is gzipped:
        filename, file_extension = os.path.splitext(readfiles[1])
        if file_extension == '.gz':
            logger.debug(f'{"[NOTE]:":10} Distributing reads from gzipped file {os.path.basename(readfiles[1])}')
            iterator2 = FastqGeneralIterator(gzip.open(readfiles[1], 'rt'))
        else:
            iterator2 = FastqGeneralIterator(open(readfiles[1]))

        for ID1_long, Seq1, Qual1 in progressbar.progressbar(iterator1, max_value=num_reads_in_readfile,
                                                             min_poll_interval=30):
            ID2_long, Seq2, Qual2 = next(iterator2)
            ID1 = ID1_long.split()[0]
            if ID1.endswith('/1') or ID1.endswith('/2'):
                ID1 = ID1[:-2]
            ID2 = ID2_long.split()[0]
            if ID2.endswith('/1') or ID2.endswith('/2'):
                ID2 = ID2[:-2]
            if ID1 in read_hit_dict:
                for target in read_hit_dict[ID1]:
                    write_paired_seqs(target, ID1, Seq1, Qual1, ID2, Seq2, Qual2, merged=merged)
                    # Note that read pairs can get written to multiple targets
            elif ID2 in read_hit_dict:
                for target in read_hit_dict[ID2]:
                    write_paired_seqs(target, ID1, Seq1, Qual1, ID2, Seq2, Qual2, merged=merged)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('bam_filename', help='Name of bam file from BWA mapping of reads')
    parser.add_argument('readfiles', nargs='+', help='List of readfiles')
    parser.add_argument("--merged", help="If provided, write fastq files for bbmerge", action="store_true",
                        default=False)
    args = parser.parse_args()

    logging.info(f'{"[NOTE]:":10} Running script distribute_reads_to_targets_bwa.py with {args}')
    readfiles = args.readfiles
    logging.info(f'{"[NOTE]:":10} readfiles are {readfiles}')
    read_hit_dict = read_sorting(args.bam_filename)
    logging.info(f'{"[NOTE]:":10} [NOTE]: Unique reads with hits: {len(read_hit_dict)}')
    distribute_reads(readfiles, read_hit_dict, merged=args.merged)


if __name__ == '__main__':
    main()
