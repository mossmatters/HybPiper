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
from collections import defaultdict
from hybpiper import utils

# Create logger:
logger = logging.getLogger(f'hybpiper.assemble.{__name__}')

# Set widget format for progressbar:
widgets = [' ' * 11,
           progressbar.Timer(),
           progressbar.Bar(),
           progressbar.ETA()]


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


def read_sorting_bwa(bamfilename):
    """
    Returns a dictionary of read_hit_dict[readID] = [target1, target2, ...]

    :param str bamfilename: path the *.bam file output by BWA
    :return: dict read_hit_dict: dictionary of read_hit_dict[readID] = [target1, target2, ...]
    """

    logger.info(f'{"[INFO]:":10} Gathering IDs for mapped reads...')

    samtools_cmd = f'samtools view -F 4 {bamfilename}'

    bwa_results = subprocess.Popen(samtools_cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)

    read_hit_dict = {}
    for line in bwa_results.stdout:
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


def read_sorting_blastx(blastfilename):
    """
    Returns a dictionary of read_hit_dict[readID] = [target1, target2, ...]

    :param str blastfilename: path the BLASTx tabular output file
    :return: dict read_hit_dict: dictionary of read_hit_dict[readID] = [target1, target2, ...]
    """

    logger.info(f'{"[INFO]:":10} Gathering IDs for mapped reads...')

    read_hit_dict = {}
    with open(blastfilename) as blastfilename_handle:
        for line in blastfilename_handle:
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
    if single:  # If True, write paired reads in interleaved format; single is historical and should be removed
        if merged:
            outfile = open(os.path.join(target, f'{target}_interleaved.fastq'), 'a')
            outfile.write(f'@{ID1}\n{Seq1}\n+\n{Qual1}\n')
            outfile.write(f'@{ID2}\n{Seq2}\n+\n{Qual2}\n')
            outfile.close()

        outfile = open(os.path.join(target, f'{target}_interleaved.fasta'), 'a')
        outfile.write(f'>{ID1}\n{Seq1}\n')
        outfile.write(f'>{ID2}\n{Seq2}\n')
        outfile.close()


def write_paired_seqs_once(target, read_list, merged=False):
    """
    Take a target name and a corresponding list of interleaved paired-end reads, and writes them to file.

    :param str target: name of the target (gene)
    :param list read_list: list of interleaved paired-end reads
    :param bool merged: If True, write fastq seqs as well as fasta
    :return:
    """

    mkdir_p(target)

    if merged:
        fastq_outfile = open(os.path.join(target, f'{target}_interleaved.fastq'), 'w')

    with open(os.path.join(target, f'{target}_interleaved.fasta'), 'w') as outfile:
        for read_pair in read_list:
            ID1 = read_pair[0]
            Seq1 = read_pair[1]
            Qual1 = read_pair[2]
            ID2 = read_pair[3]
            Seq2 = read_pair[4]
            Qual2 = read_pair[5]

            outfile.write(f'>{ID1}\n{Seq1}\n')
            outfile.write(f'>{ID2}\n{Seq2}\n')

            if merged:
                fastq_outfile.write(f'>{ID1}\n{Seq1}\n+\n{Qual1}\n')
                fastq_outfile.write(f'>{ID2}\n{Seq2}\n+\n{Qual2}\n')

    if merged:
        fastq_outfile.close()


def write_single_seqs(target, ID1, Seq1):
    """
    Writes a fasta file of single-end/unpaired reads to the corresponding gene directory

    :param str target: gene name e.g. gene001
    :param str ID1: fasta/fastq header for R1
    :param str Seq1: fasta/fastq sequence for R1
    :return::
    """

    mkdir_p(target)
    outfile = open(os.path.join(target, f'{target}_unpaired.fasta'), 'a')
    outfile.write(f'>{ID1}\n{Seq1}\n')
    outfile.close()


def write_single_seqs_once(target, read_list):
    """
    Take a target name and a corresponding list of single-end reads, and writes them to file.

    :param str target: name of the target (gene)
    :param list read_list: list of interleaved paired-end reads
    :return:
    """

    mkdir_p(target)

    with open(os.path.join(target, f'{target}_unpaired.fasta'), 'w') as outfile:
        for read_pair in read_list:
            ID1 = read_pair[0]
            Seq1 = read_pair[1]
            outfile.write(f'>{ID1}\n{Seq1}\n')


def distribute_reads(readfiles, read_hit_dict, merged=False, unpaired_readfile=None, single_end=False, low_mem=False):
    """

    :param list readfiles: a list of one or more readfiles
    :param dict read_hit_dict: dictionary of read_hit_dict[readID] = [target1, target2, ...]
    :param bool merged: boolean passed to function write_paired_seqs()
    :param str/bool unpaired_readfile: a path if an unpaired file has been provided, False if not
    :param bool single_end: True if a single file was provided as input to -r, False if not
    :param bool low_mem: If False, reads to distribute will be saved in a dictionary and written once; uses more RAM
    :return:
    """

    if merged:
        logger.info(f'{"[NOTE]:":10} Writing fastq files for merging with BBmerge.sh')

    if low_mem:
        logger.info(f'{"[NOTE]:":10} Read distribution running in low-mem mode; note that this can use less '
                    f'memory/RAM, but is slower!')

    gene_2_reads_dict = defaultdict(list)  # created regardless of low_mem mode

    num_reads_in_readfile = 0

    # Check if read file is gzipped:
    filename, file_extension = os.path.splitext(readfiles[0])
    if file_extension == '.gz':
        logger.debug(f'Distributing reads from gzipped file {os.path.basename(readfiles[0])}')
        for line in gzip.open(readfiles[0], 'rt'):
            num_reads_in_readfile += 1  # Get total # reads for progressbar and to write file for hybpiper_stats.py
        num_reads_in_readfile = int(num_reads_in_readfile / 4)
        iterator1 = FastqGeneralIterator(gzip.open(readfiles[0], 'rt'))

    else:
        for line in open(readfiles[0]):
            num_reads_in_readfile += 1  # Get total # reads for progressbar and to write file for hybpiper_stats.py
        num_reads_in_readfile = int(num_reads_in_readfile / 4)
        iterator1 = FastqGeneralIterator(open(readfiles[0]))

    if len(readfiles) == 1 and single_end:
        logger.info(f'{"[NOTE]:":10} Distributing single-end reads to gene directories')

        try:
            for ID1_long, Seq1, Qual1 in progressbar.progressbar(iterator1, max_value=num_reads_in_readfile,
                                                                 min_poll_interval=30, widgets=widgets):
                ID1 = ID1_long.split()[0]
                if ID1.endswith('\1') or ID1.endswith('\2'):
                    ID1 = ID1[:-2]
                if ID1 in read_hit_dict:
                    for target in read_hit_dict[ID1]:
                        if low_mem:
                            write_single_seqs(target, ID1, Seq1)
                        else:
                            gene_2_reads_dict[target].append((ID1, Seq1))
        except ValueError as e:
            fill = utils.fill_forward_slash(f'')
            logger.error(fill)
            fill = utils.fill_forward_slash(f'{"[ERROR]:":10} Malformed FASTQ file! Please check your input FASTQ file '
                                            f'({readfiles[0]}) and try again. Error is: {e}',
                                            width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                            break_on_forward_slash=True)
            logger.error(f'{fill}')
            sys.exit()

        if not low_mem:
            for target, read_list in gene_2_reads_dict.items():
                write_single_seqs_once(target, read_list)

        # Write a file containing the total number of single-end reads in the input file, to be parsed by
        # hybpiper_stats.py when calculating BLASTX enrichment efficiency:
        with open(f'total_input_reads_single.txt', 'w') as single_reads_number:
            single_reads_number.write(f'{num_reads_in_readfile}\n')

    if len(readfiles) == 1 and unpaired_readfile:
        logger.info(f'{"[NOTE]:":10} Distributing unpaired reads to gene directories')

        try:
            for ID1_long, Seq1, Qual1 in progressbar.progressbar(iterator1, max_value=num_reads_in_readfile,
                                                                 min_poll_interval=30, widgets=widgets):
                ID1 = ID1_long.split()[0]
                if ID1.endswith('\1') or ID1.endswith('\2'):
                    ID1 = ID1[:-2]
                if ID1 in read_hit_dict:
                    for target in read_hit_dict[ID1]:
                        if low_mem:
                            write_single_seqs(target, ID1, Seq1)
                        else:
                            gene_2_reads_dict[target].append((ID1, Seq1))
        except ValueError as e:
            fill = utils.fill_forward_slash(f'')
            logger.error(fill)
            fill = utils.fill_forward_slash(f'{"[ERROR]:":10} Malformed FASTQ file! Please check your input FASTQ file '
                                            f'({readfiles[0]}) and try again. Error is: {e}',
                                            width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                            break_on_forward_slash=True)
            logger.error(f'{fill}')
            sys.exit()

        if not low_mem:
            for target, read_list in gene_2_reads_dict.items():
                write_single_seqs_once(target, read_list)

        # Write a file containing the total number of unpaired reads in the input file, to be parsed by
        # hybpiper_stats.py when calculating BLASTX enrichment efficiency:
        with open(f'total_input_reads_unpaired.txt', 'w') as unpaired_reads_number:
            unpaired_reads_number.write(f'{num_reads_in_readfile}\n')

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

        try:
            for ID1_long, Seq1, Qual1 in progressbar.progressbar(iterator1, max_value=num_reads_in_readfile,
                                                                 min_poll_interval=30, widgets=widgets):
                ID2_long, Seq2, Qual2 = next(iterator2)
                ID1 = ID1_long.split()[0]
                if ID1.endswith('/1') or ID1.endswith('/2'):
                    ID1 = ID1[:-2]
                ID2 = ID2_long.split()[0]
                if ID2.endswith('/1') or ID2.endswith('/2'):
                    ID2 = ID2[:-2]
                if ID1 in read_hit_dict:
                    for target in read_hit_dict[ID1]:
                        if low_mem:
                            write_paired_seqs(target, ID1, Seq1, Qual1, ID2, Seq2, Qual2, merged=merged)
                            # Note that read pairs can get written to multiple targets
                        else:
                            gene_2_reads_dict[target].append((ID1, Seq1, Qual1, ID2, Seq2, Qual2))
                elif ID2 in read_hit_dict:
                    for target in read_hit_dict[ID2]:
                        if low_mem:
                            write_paired_seqs(target, ID1, Seq1, Qual1, ID2, Seq2, Qual2, merged=merged)
                        else:
                            gene_2_reads_dict[target].append((ID1, Seq1, Qual1, ID2, Seq2, Qual2))
        except ValueError as e:
            fill = utils.fill_forward_slash(f'')
            logger.error(fill)
            fill = utils.fill_forward_slash(f'{"[ERROR]:":10} Malformed FASTQ file! Please check your input FASTQ '
                                            f'files ({readfiles[0]}, {readfiles[1]}) and try again. Error is: {e}',
                                            width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                            break_on_forward_slash=True)
            logger.error(f'{fill}')
            sys.exit()

        if not low_mem:
            for target, read_list in gene_2_reads_dict.items():
                write_paired_seqs_once(target, read_list, merged=merged)

        # Write a file containing the total number of unpaired reads in the input file, to be parsed by
        # hybpiper_stats.py when calculating BLASTX enrichment efficiency:
        with open(f'total_input_reads_paired.txt', 'w') as paired_reads_number:
            paired_reads_number.write(f'{num_reads_in_readfile * 2}\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('bam_filename', help='Name of bam file from BWA mapping of reads', default=False)
    group.add_argument('blast_filename', help='Name of blast file from BLASTX mapping of reads', default=False)
    parser.add_argument('readfiles', nargs='+', help='List of readfiles')
    parser.add_argument("--merged", help="If provided, write fastq files for bbmerge", action="store_true",
                        default=False)
    args = parser.parse_args()

    logging.info(f'{"[NOTE]:":10} Running script distribute_reads_to_targets_bwa.py with {args}')
    readfiles = args.readfiles
    logging.info(f'{"[NOTE]:":10} readfiles are {readfiles}')

    if args.bam_filename:
        read_hit_dict = read_sorting_bwa(args.bam_filename)
    else:
        read_hit_dict = read_sorting_blastx(args.blast_filename)

    logging.info(f'{"[NOTE]:":10} [NOTE]: Unique reads with hits: {len(read_hit_dict)}')
    distribute_reads(readfiles, read_hit_dict, merged=args.merged)


if __name__ == '__main__':
    main()
