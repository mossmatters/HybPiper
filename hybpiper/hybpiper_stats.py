#!/usr/bin/env python

#########################
# HybPiper stats script #
#########################

"""
Writes a report file called "seq_lengths.tsv" (default filename, user can change this). The first line contains gene
names. The second line contains the length of the reference sequences (targets). If there are multiple targets per gene,
the mean length is reported. All other rows contain one sample per line.

Parses the "seq_lengths.tsv" and gathers additional statistics about the HybPiper run. Writes a report file called
"hybpiper_stats.tsv" (default filename, user can change this).

For an explanation of columns, see github.com/mossmatters/HybPiper/wiki
"""

import argparse
import os
import sys
import subprocess
import re
from Bio import SeqIO
from collections import defaultdict
import logging
import textwrap
import progressbar
import multiprocessing
from multiprocessing import Manager
from concurrent.futures.process import ProcessPoolExecutor
from concurrent.futures import wait, as_completed
import traceback
import tarfile
import io

from hybpiper import utils
from hybpiper.version import __version__


# Create a custom logger

# Log to Terminal (stderr):
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)

# Setup logger:
logger = logging.getLogger(f'hybpiper.{__name__}')

# Add handlers to the logger
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)  # Default level is 'WARNING'

# Set widget format for progressbar:
widgets = [
    ' ' * 11,
    progressbar.Timer(),
    progressbar.Bar(),
    progressbar.ETA()
]


def seq_length_calc(seq_lengths_fn):
    """
    From the output file produced by get_seq_lengths(), calculate the number of genes with seqs, and at least a
    percentage of the reference length.

    :param path seq_lengths_fn: path to the sequence length file produced by get_seq_lengths()
    :return dict seq_length_dict: dictionary of sample:list_of_genes_above_length_threshold
    """

    seq_length_dict = {}

    with open(seq_lengths_fn) as seq_len:
        gene_names = seq_len.readline()  # skip the first line
        target_lengths = seq_len.readline().split()[1:]

        for line in seq_len:
            line = line.split()
            name = line.pop(0)
            is_25pct = 0
            is_50pct = 0
            is_75pct = 0
            is_150pct = 0

            for gene in range(len(line)):
                gene_length = float(line[gene])
                target_length = float(target_lengths[gene])
                if gene_length > target_length * 0.25:
                    is_25pct += 1
                if gene_length > target_length * 0.50:
                    is_50pct += 1
                if gene_length > target_length * 0.75:
                    is_75pct += 1
                if gene_length > target_length * 1.5:
                    is_150pct += 1

            seq_length_dict[name] = [str(is_25pct), str(is_50pct), str(is_75pct), str(is_150pct)]

    return seq_length_dict


def parse_sample(sample_name,
                 sampledir_parent,
                 compressed_samples_set,
                 compressed_sample_dict,
                 unique_gene_names,
                 sequence_type,
                 lock,
                 counter):
    """

    :param str sample_name: name of the sample (no ".tar.gz")
    :param str sampledir_parent: path of the parent directory containing HybPiper output
    :param set compressed_samples_set: set of sample names that are compressed (i.e. *.tar.gz)
    :param dict compressed_sample_dict: dict of compressed file contents for compressed samples
    :param list unique_gene_names: a list of unique gene/locus names extracted from the target file
    :param str sequence_type: sequence type (gene or supercontig) to recover stats for
    :param multiprocessing.managers.AcquirerProxy lock:
    :param multiprocessing.managers.ValueProxy counter:
    :return:
    """

    compressed_bool = True if sample_name in compressed_samples_set else False
    warning_messages_list = []

    ####################################################################################################################
    # Set paths of files to get/process for the sample:
    ####################################################################################################################

    gene_fasta_paths = []
    if sequence_type.upper() == 'GENE':
        for gene in unique_gene_names:
            seq_file = os.path.join(sample_name,
                                    gene,
                                    sample_name,
                                    'sequences',
                                    'FNA',
                                    f'{gene}.FNA')
            gene_fasta_paths.append(seq_file)

    else:
        for gene in unique_gene_names:
            seq_file = os.path.join(sample_name,
                                    gene,
                                    sample_name,
                                    'sequences',
                                    'intron',
                                    f'{gene}_supercontig.FNA')
            gene_fasta_paths.append(seq_file)

    bamfile_fn = f'{sample_name}/{sample_name}.bam'
    bamfile_unpaired_fn = f'{sample_name}/{sample_name}_unpaired.bam'
    bam_flagstats_tsv_fn = f'{sample_name}/{sample_name}_bam_flagstat.tsv'
    unpaired_bam_flagstats_tsv_fn = f'{sample_name}/{sample_name}_unpaired_bam_flagstat.tsv'

    blastxfile_fn = f'{sample_name}/{sample_name}.blastx'
    blastxfile_unpaired_fn = f'{sample_name}/{sample_name}_unpaired.blastx'

    total_input_reads_paired_fn = f'{sample_name}/total_input_reads_paired.txt'
    total_input_reads_single_fn = f'{sample_name}/total_input_reads_single.txt'
    total_input_reads_unpaired_fn = f'{sample_name}/total_input_reads_unpaired.txt'

    loci_with_mapped_reads_fn = f'{sample_name}/spades_genelist.txt'
    loci_with_contigs_fn = f'{sample_name}/exonerate_genelist.txt'
    loci_with_extracted_seqs_fn = f'{sample_name}/genes_with_seqs.txt'

    long_paralog_warnings_file_fn = f'{sample_name}/{sample_name}_genes_with_long_paralog_warnings.txt'
    depth_paralog_warnings_file_fn = f'{sample_name}/{sample_name}_genes_with_paralog_warnings_by_contig_depth.csv'

    genes_with_stitched_contig_file_fn = f'{sample_name}/{sample_name}_genes_with_stitched_contig.csv'

    genes_derived_from_putative_chimeric_stitched_contig_file_fn = \
        f'{sample_name}/{sample_name}_genes_derived_from_putative_chimeric_stitched_contig.csv'

    ####################################################################################################################
    # Check that a bam or blastx mapping file is present, and check if `samtools stats` has already been run if bam:
    ####################################################################################################################
    no_mapping_file = False
    if compressed_bool:
        if (bamfile_fn not in compressed_sample_dict[sample_name] and
                blastxfile_fn not in compressed_sample_dict[sample_name]):
            no_mapping_file = True
    elif (not utils.file_exists_and_not_empty(f'{sampledir_parent}/{bamfile_fn}') and
          not utils.file_exists_and_not_empty(f'{sampledir_parent}/{blastxfile_fn}')):
        no_mapping_file = True

    if no_mapping_file:
        warning_messages_list.append(f'No *.bam or *.blastx file found for {sample_name}. No statistics will be '
                                     f'recovered. Please check the log file for this sample!')

        with lock:
            counter.value += 1
            return (sample_name,
                    None,
                    counter.value,
                    warning_messages_list)

    bam_flagstats_tsv_file_bool = False
    if compressed_bool:
        if (bamfile_fn in compressed_sample_dict[sample_name] and bam_flagstats_tsv_fn in
                compressed_sample_dict[sample_name]):
            bam_flagstats_tsv_file_bool = True
    else:
        if (utils.file_exists_and_not_empty(f'{sampledir_parent}/{bamfile_fn}') and
                utils.file_exists_and_not_empty(f'{sampledir_parent}/{bam_flagstats_tsv_fn}')):
            bam_flagstats_tsv_file_bool = True

    ####################################################################################################################
    # Check that a file containing total input paired or single-end read count can be found:
    ####################################################################################################################
    no_count_file = False
    if compressed_bool:
        if (total_input_reads_paired_fn not in compressed_sample_dict[sample_name] and
                total_input_reads_single_fn not in compressed_sample_dict[sample_name]):
            no_count_file = True
    elif not (utils.file_exists_and_not_empty(f'{sampledir_parent}/{total_input_reads_paired_fn}') and
              not utils.file_exists_and_not_empty(f'{sampledir_parent}/{total_input_reads_single_fn}')):
        no_count_file = True

    if no_count_file:
        warning_messages_list.append(f'No file containing total input paired or single-end read count found for '
                                     f'sample {sample_name}. No statistics will be recovered. Please check the log '
                                     f'file for this sample!')

        with lock:
            counter.value += 1
            return (sample_name,
                    None,
                    counter.value,
                    warning_messages_list)

    ####################################################################################################################
    # Set some default stats:
    ####################################################################################################################
    sample_stats_dict = dict()

    sample_stats_dict['seq_lengths'] = dict()
    for gene_name in unique_gene_names:
        sample_stats_dict['seq_lengths'][gene_name] = 0  # set default of zero length

    sample_stats_dict['sample_total_bases'] = 0

    total_input_reads_paired = 0
    total_input_reads_single = 0
    total_input_reads_unpaired = 0

    mapped_reads_main = 0
    mapped_reads_unpaired = 0

    loci_with_mapped_reads = 0
    loci_with_contigs = 0
    loci_with_extracted_seqs = 0

    paralog_warnings_long = 0
    num_genes_paralog_warning_by_depth = 0

    stitched_contig_produced = 0
    no_stitched_contig = 0
    stitched_contig_skipped = 0

    chimeric_stitched_contigs = 0

    ####################################################################################################################
    # Iterate over the files/folders for the sample and recover stats:
    ####################################################################################################################
    if compressed_bool:
        compressed_sample_file_path = f'{sampledir_parent}/{sample_name}.tar.gz'

        with tarfile.open(compressed_sample_file_path, 'r:gz') as tarfile_handle:

            for tarinfo in tarfile_handle.getmembers():
                if tarinfo.size != 0:

                    ####################################################################################################
                    # Get seq lengths:
                    ####################################################################################################
                    if tarinfo.name in gene_fasta_paths:
                        gene_name, _ = os.path.splitext(os.path.basename(tarinfo.name))

                        seqrecord = utils.get_compressed_seqrecords(tarfile_handle,
                                                                    tarinfo)[0]
                        seq_length = len(seqrecord.seq.replace('N', ''))
                        sample_stats_dict['seq_lengths'][gene_name] = seq_length
                        sample_stats_dict['sample_total_bases'] += seq_length

                    ####################################################################################################
                    # Get bam mapping stats if present, or blastx if not:
                    ####################################################################################################
                    if tarinfo.name == bamfile_fn and not bam_flagstats_tsv_file_bool:
                        # Run samtools flagstats now
                        bam_flagstats_lines, message_list = \
                            utils.get_bamtools_flagstat_lines_from_compressed(sample_name,
                                                                              sampledir_parent,
                                                                              bamfile_fn)
                        for line in bam_flagstats_lines:
                            if re.search(r'\bprimary mapped$\b', line):
                                mapped_reads_main += int(line.split('\t')[0])

                        for message in message_list:
                            warning_messages_list.append(message)

                    if tarinfo.name == bamfile_unpaired_fn and not bam_flagstats_tsv_file_bool:
                        # Run samtools flagstats now
                        unpaired_bam_flagstats_lines, message_list = \
                            utils.get_bamtools_flagstat_lines_from_compressed(sample_name,
                                                                              sampledir_parent,
                                                                              bamfile_unpaired_fn)
                        for line in unpaired_bam_flagstats_lines:
                            if re.search(r'\bprimary mapped$\b', line):
                                mapped_reads_unpaired += int(line.split('\t')[0])

                        for message in message_list:
                            warning_messages_list.append(message)

                    if tarinfo.name == bam_flagstats_tsv_fn:
                        bam_flagstats_lines = utils.get_compressed_file_lines(tarfile_handle,
                                                                              tarinfo)

                        for line in bam_flagstats_lines:
                            if re.search(r'\bprimary mapped$\b', line):
                                mapped_reads_main += int(line.split('\t')[0])

                    if tarinfo.name == unpaired_bam_flagstats_tsv_fn:
                        unpaired_bam_flagstats_lines = utils.get_compressed_file_lines(tarfile_handle,
                                                                                       tarinfo)
                        for line in unpaired_bam_flagstats_lines:
                            if re.search(r'\bprimary mapped$\b', line):
                                mapped_reads_unpaired += int(line.split('\t')[0])

                    if tarinfo.name == blastxfile_fn and bamfile_fn not in compressed_sample_dict[sample_name]:
                        blastx_lines_main = utils.get_compressed_file_lines(tarfile_handle,
                                                                            tarinfo)
                        reads_with_hits_main = [x.split('\t')[0] for x in blastx_lines_main if x]
                        mapped_reads_main = len(set(reads_with_hits_main))

                    if tarinfo.name == blastxfile_unpaired_fn and bamfile_fn not in compressed_sample_dict[sample_name]:
                        blastx_lines_unpaired = utils.get_compressed_file_lines(tarfile_handle,
                                                                                tarinfo)
                        reads_with_hits_unpaired = [x.split('\t')[0] for x in blastx_lines_unpaired if x]
                        mapped_reads_unpaired = len(set(reads_with_hits_unpaired))

                    ####################################################################################################
                    # Get input read stats
                    ####################################################################################################
                    if tarinfo.name == total_input_reads_paired_fn:
                        total_input_reads_paired_lines = utils.get_compressed_file_lines(tarfile_handle,
                                                                                         tarinfo)
                        total_input_reads_paired = int(total_input_reads_paired_lines[0].rstrip())

                    if tarinfo.name == total_input_reads_single_fn:
                        total_input_reads_single_lines = utils.get_compressed_file_lines(tarfile_handle,
                                                                                         tarinfo)
                        total_input_reads_single = int(total_input_reads_single_lines[0].rstrip())

                    if tarinfo.name == total_input_reads_unpaired_fn:
                        total_input_reads_unpaired_lines = utils.get_compressed_file_lines(tarfile_handle,
                                                                                           tarinfo)
                        total_input_reads_unpaired = int(total_input_reads_unpaired_lines[0].rstrip())

                    ####################################################################################################
                    # Get loci with mapped reads:
                    ####################################################################################################
                    if tarinfo.name == loci_with_mapped_reads_fn:
                        loci_with_mapped_reads_lines = utils.get_compressed_file_lines(tarfile_handle,
                                                                                       tarinfo)
                        loci_with_mapped_reads = len(loci_with_mapped_reads_lines)

                    ####################################################################################################
                    # Get loci with contigs:
                    ####################################################################################################
                    if tarinfo.name == loci_with_contigs_fn:
                        loci_with_contigs_lines = utils.get_compressed_file_lines(tarfile_handle,
                                                                                  tarinfo)
                        loci_with_contigs = len(loci_with_contigs_lines)

                    ####################################################################################################
                    # Get loci with extracted seqs:
                    ####################################################################################################
                    if tarinfo.name == loci_with_extracted_seqs_fn:
                        loci_with_extracted_seqs_lines = utils.get_compressed_file_lines(tarfile_handle,
                                                                                         tarinfo)
                        loci_with_extracted_seqs = len(loci_with_extracted_seqs_lines)

                    ####################################################################################################
                    # Get paralog warnings - long:
                    ####################################################################################################
                    if tarinfo.name == long_paralog_warnings_file_fn:
                        long_paralog_warnings_file_lines = utils.get_compressed_file_lines(tarfile_handle,
                                                                                           tarinfo)
                        paralog_warnings_long = len(long_paralog_warnings_file_lines)

                    ####################################################################################################
                    # Get paralogs warnings - by contig depth:
                    ####################################################################################################
                    if tarinfo.name == depth_paralog_warnings_file_fn:
                        depth_paralog_warnings_file_lines = utils.get_compressed_file_lines(tarfile_handle,
                                                                                            tarinfo)
                        for gene_stats in depth_paralog_warnings_file_lines:
                            stat = gene_stats.split(',')[3].strip()
                            if stat == 'True':
                                num_genes_paralog_warning_by_depth += 1

                    ####################################################################################################
                    # Stitched contig stats:
                    ####################################################################################################
                    if tarinfo.name == genes_with_stitched_contig_file_fn:
                        stitched_contig_lines = utils.get_compressed_file_lines(tarfile_handle,
                                                                                tarinfo)
                        for gene_stats in stitched_contig_lines:
                            stat = gene_stats.split(',')[2]
                            if re.search('single Exonerate hit', stat):
                                no_stitched_contig += 1
                            elif re.search('Stitched contig produced', stat):
                                stitched_contig_produced += 1
                            elif re.search('Stitched contig step skipped', stat):
                                stitched_contig_skipped += 1

                    ####################################################################################################
                    # Chimeric contigs stats:
                    ####################################################################################################
                    if tarinfo.name == genes_derived_from_putative_chimeric_stitched_contig_file_fn:

                        chimeric_contigs_lines = utils.get_compressed_file_lines(tarfile_handle,
                                                                                 tarinfo)

                        for gene_stats in chimeric_contigs_lines:
                            stat = gene_stats.split(',')[2]
                            if re.search(' Chimera WARNING for stitched_contig.', stat):
                                chimeric_stitched_contigs += 1

    else:  # i.e. sample is uncompressed

        ################################################################################################################
        # Get seq lengths:
        ################################################################################################################
        for gene_path in gene_fasta_paths:
            gene_name, _ = os.path.splitext(os.path.basename(gene_path))
            gene_path_full = f'{sampledir_parent}/{gene_path}'

            if utils.file_exists_and_not_empty(gene_path_full):
                seqrecord = SeqIO.read(gene_path_full, 'fasta')
                seq_length = len(seqrecord.seq.replace('N', ''))
                sample_stats_dict['seq_lengths'][gene_name] = seq_length
                sample_stats_dict['sample_total_bases'] += seq_length

        ################################################################################################################
        # Get bam mapping stats if present, or blastx if not:
        ################################################################################################################
        bamfile_fn_full_path = f'{sampledir_parent}/{bamfile_fn}'
        bamfile_unpaired_fn_full_path = f'{sampledir_parent}/{bamfile_unpaired_fn}'
        blastx_fn_full_path = f'{sampledir_parent}/{blastxfile_fn}'
        blastxfile_unpaired_fn_full_path = f'{sampledir_parent}/{blastxfile_unpaired_fn}'
        bam_flagstats_tsv_fn_full_path = f'{sampledir_parent}/{bam_flagstats_tsv_fn}'
        unpaired_bam_flagstats_tsv_fn_full_path = f'{sampledir_parent}/{unpaired_bam_flagstats_tsv_fn}'

        if utils.file_exists_and_not_empty(bamfile_fn_full_path) and not bam_flagstats_tsv_file_bool:
            # Run samtools flagstats now:
            bam_flagstats_lines, message_list = \
                utils.get_bamtools_flagstat_lines_from_uncompressed(sample_name,
                                                                    sampledir_parent,
                                                                    bamfile_fn)
            for line in bam_flagstats_lines:
                if re.search(r'\bprimary mapped$\b', line):
                    mapped_reads_main += int(line.split('\t')[0])

            for message in message_list:
                warning_messages_list.append(message)

        if utils.file_exists_and_not_empty(bamfile_unpaired_fn_full_path) and not bam_flagstats_tsv_file_bool:
            # Run samtools flagstats now:
            unpaired_bam_flagstats_lines, message_list = \
                utils.get_bamtools_flagstat_lines_from_uncompressed(sample_name,
                                                                    sampledir_parent,
                                                                    bamfile_unpaired_fn)
            for line in unpaired_bam_flagstats_lines:
                if re.search(r'\bprimary mapped$\b', line):
                    mapped_reads_unpaired += int(line.split('\t')[0])

            for message in message_list:
                warning_messages_list.append(message)

        if utils.file_exists_and_not_empty(bam_flagstats_tsv_fn_full_path):
            with open(bam_flagstats_tsv_fn_full_path) as bam_flagstats_tsv_fn_full_path_handle:
                bam_flagstats_lines = list(bam_flagstats_tsv_fn_full_path_handle.readlines())
                for line in bam_flagstats_lines:
                    if re.search(r'\bprimary mapped$\b', line):
                        mapped_reads_main += int(line.split('\t')[0])

        if utils.file_exists_and_not_empty(unpaired_bam_flagstats_tsv_fn_full_path):
            with open(unpaired_bam_flagstats_tsv_fn_full_path) as unpaired_bam_flagstats_tsv_fn_full_path_handle:
                unpaired_bam_flagstats_lines = list(unpaired_bam_flagstats_tsv_fn_full_path_handle.readlines())
                for line in unpaired_bam_flagstats_lines:
                    if re.search(r'\bprimary mapped$\b', line):
                        mapped_reads_unpaired += int(line.split('\t')[0])

        if (utils.file_exists_and_not_empty(blastx_fn_full_path) and
                not utils.file_exists_and_not_empty(bamfile_fn_full_path)):

            with open(blastx_fn_full_path) as blastx_fn_handle:
                blastx_lines_main = blastx_fn_handle.readlines()
                reads_with_hits_main = [x.split('\t')[0] for x in blastx_lines_main if x]
                mapped_reads_main = len(set(reads_with_hits_main))

        if (utils.file_exists_and_not_empty(blastxfile_unpaired_fn_full_path) and
                not utils.file_exists_and_not_empty(bamfile_fn_full_path)):

            with open(blastxfile_unpaired_fn_full_path) as blastxfile_unpaired_fn_handle:
                blastx_lines_unpaired = blastxfile_unpaired_fn_handle.readlines()
                reads_with_hits_unpaired = [x.split('\t')[0] for x in blastx_lines_unpaired if x]
                mapped_reads_unpaired = len(set(reads_with_hits_unpaired))

        ################################################################################################################
        # Get input read stats
        ################################################################################################################
        total_input_reads_paired_fn_full_path = f'{sampledir_parent}/{total_input_reads_paired_fn}'
        if utils.file_exists_and_not_empty(total_input_reads_paired_fn_full_path):
            with open(total_input_reads_paired_fn_full_path) as paired_number_handle:
                total_input_reads_paired = int(paired_number_handle.read().rstrip())

        total_input_reads_single_fn_full_path = f'{sampledir_parent}/{total_input_reads_single_fn}'
        if utils.file_exists_and_not_empty(total_input_reads_single_fn_full_path):
            with open(total_input_reads_single_fn_full_path) as single_number_handle:
                total_input_reads_single = int(single_number_handle.read().rstrip())

        total_input_reads_unpaired_fn_full_path = f'{sampledir_parent}/{total_input_reads_unpaired_fn}'
        if utils.file_exists_and_not_empty(total_input_reads_unpaired_fn_full_path):
            with open(total_input_reads_unpaired_fn_full_path) as unpaired_number_handle:
                total_input_reads_unpaired = int(unpaired_number_handle.read().rstrip())

        ################################################################################################################
        # Get loci with mapped reads:
        ################################################################################################################
        loci_with_mapped_reads_fn_full_path = f'{sampledir_parent}/{loci_with_mapped_reads_fn}'
        if utils.file_exists_and_not_empty(loci_with_mapped_reads_fn_full_path):
            with open(loci_with_mapped_reads_fn_full_path) as mapped_reads_handle:
                loci_with_mapped_reads_lines = mapped_reads_handle.readlines()
                loci_with_mapped_reads = len(loci_with_mapped_reads_lines)

        ################################################################################################################
        # Get loci with contigs:
        ################################################################################################################
        loci_with_contigs_fn_full_path = f'{sampledir_parent}/{loci_with_contigs_fn}'
        if utils.file_exists_and_not_empty(loci_with_contigs_fn_full_path):
            with open(loci_with_contigs_fn_full_path) as contigs_handle:
                loci_with_contigs_lines = contigs_handle.readlines()
                loci_with_contigs = len(loci_with_contigs_lines)

        ################################################################################################################
        # Get loci with extracted seqs:
        ################################################################################################################
        loci_with_extracted_seqs_fn_full_path = f'{sampledir_parent}/{loci_with_extracted_seqs_fn}'
        if utils.file_exists_and_not_empty(loci_with_extracted_seqs_fn_full_path):
            with open(loci_with_extracted_seqs_fn_full_path) as extracted_seqs_handle:
                loci_with_extracted_seqs_lines = extracted_seqs_handle.readlines()
                loci_with_extracted_seqs = len(loci_with_extracted_seqs_lines)

        ################################################################################################################
        # Get paralog warnings - long:
        ################################################################################################################
        long_paralog_warnings_file_fn_full_path = f'{sampledir_parent}/{long_paralog_warnings_file_fn}'
        if utils.file_exists_and_not_empty(long_paralog_warnings_file_fn_full_path):
            with open(long_paralog_warnings_file_fn_full_path) as long_paralogs_handle:
                paralog_warnings_long = len([line for line in long_paralogs_handle.readlines() if line.rstrip()])

        ################################################################################################################
        # Get paralogs warnings - by contig depth:
        ################################################################################################################
        depth_paralog_warnings_file_fn_full_path = f'{sampledir_parent}/{depth_paralog_warnings_file_fn}'
        if utils.file_exists_and_not_empty(depth_paralog_warnings_file_fn_full_path):
            with open(depth_paralog_warnings_file_fn_full_path) as depth_paralogs_handle:
                lines = depth_paralogs_handle.readlines()
                for gene_stats in lines:
                    stat = gene_stats.split(',')[3].strip()
                    if stat == 'True':
                        num_genes_paralog_warning_by_depth += 1

        ################################################################################################################
        # Stitched contig stats:
        ################################################################################################################
        depth_paralog_warnings_file_fn_full_path = f'{sampledir_parent}/{genes_with_stitched_contig_file_fn}'

        if utils.file_exists_and_not_empty(depth_paralog_warnings_file_fn_full_path):
            with open(depth_paralog_warnings_file_fn_full_path) as stitched_contig_stats_handle:
                lines = stitched_contig_stats_handle.readlines()
                for gene_stats in lines:
                    stat = gene_stats.split(',')[2]
                    if re.search('single Exonerate hit', stat):
                        no_stitched_contig += 1
                    elif re.search('Stitched contig produced', stat):
                        stitched_contig_produced += 1
                    elif re.search('Stitched contig step skipped', stat):
                        stitched_contig_skipped += 1

        ################################################################################################################
        # Chimeric contigs stats:
        ################################################################################################################
        genes_derived_from_putative_chimeric_stitched_contig_file_fn_full_path = \
            f'{sampledir_parent}/{genes_derived_from_putative_chimeric_stitched_contig_file_fn}'

        if utils.file_exists_and_not_empty(genes_derived_from_putative_chimeric_stitched_contig_file_fn_full_path):
            with (open(genes_derived_from_putative_chimeric_stitched_contig_file_fn_full_path) as
                  chimeric_stitched_contig_stats_handle):

                lines = chimeric_stitched_contig_stats_handle.readlines()
                for gene_stats in lines:
                    stat = gene_stats.split(',')[2]
                    if re.search(' Chimera WARNING for stitched_contig.', stat):
                        chimeric_stitched_contigs += 1

    ####################################################################################################################
    # Populate sample stat dict:
    ####################################################################################################################

    ####################################################################################################################
    # Calculate enrichment efficiency:
    ####################################################################################################################
    total_input_reads = total_input_reads_paired + total_input_reads_single + total_input_reads_unpaired
    total_mapped_reads = mapped_reads_main + mapped_reads_unpaired

    try:
        pct_mapped = 100 * total_mapped_reads / total_input_reads
    except ZeroDivisionError:
        pct_mapped = 0.0

    sample_stats_dict['total_input_reads'] = str(total_input_reads)
    sample_stats_dict['total_mapped_reads'] = str(total_mapped_reads)
    sample_stats_dict['pct_mapped'] = "{0:.1f}".format(pct_mapped)

    ####################################################################################################################
    # Calculate recovery efficiency:
    ####################################################################################################################
    sample_stats_dict['loci_with_mapped_reads'] = str(loci_with_mapped_reads)
    sample_stats_dict['loci_with_contigs'] = str(loci_with_contigs)
    sample_stats_dict['loci_with_extracted_seqs'] = str(loci_with_extracted_seqs)

    ####################################################################################################################
    # Paralogs - long:
    ####################################################################################################################
    sample_stats_dict['paralog_warnings_long'] = str(paralog_warnings_long)

    ####################################################################################################################
    # Paralogs - by contig depth across query protein:
    ####################################################################################################################
    sample_stats_dict['paralog_warnings_depth'] = str(num_genes_paralog_warning_by_depth)

    ####################################################################################################################
    # Stitched contig information:
    ####################################################################################################################
    sample_stats_dict['no_stitched_contig'] = str(no_stitched_contig)
    sample_stats_dict['stitched_contig_produced'] = str(stitched_contig_produced)
    sample_stats_dict['stitched_contig_skipped'] = str(stitched_contig_skipped)

    ####################################################################################################################
    # Chimeric stitched contigs:
    ####################################################################################################################
    sample_stats_dict['chimeric_contig_warnings'] = str(chimeric_stitched_contigs)

    with (lock):
        counter.value += 1

        return (sample_name,
                sample_stats_dict,
                counter.value,
                warning_messages_list)


def standalone():
    """
    Used when this module is run as a stand-alone script. Parses command line arguments and runs function main().:
    """

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    group_1 = parser.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna',
                         dest='targetfile_dna',
                         default=False,
                         help='FASTA file containing DNA target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa',
                         dest='targetfile_aa',
                         default=False,
                         help='FASTA file containing amino-acid target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    parser.add_argument("sequence_type",
                        choices=["gene", "GENE", "supercontig", "SUPERCONTIG"],
                        help="Sequence type (gene or supercontig) to recover stats for.")
    parser.add_argument("namelist",
                        help="text file with names of HybPiper output directories, one per line.")
    parser.add_argument("--seq_lengths_filename",
                        help="File name for the sequence lengths *.tsv file. Default is <seq_lengths.tsv>.",
                        default='seq_lengths')
    parser.add_argument("--stats_filename",
                        help="File name for the stats *.tsv file. Default is <hybpiper_stats.tsv>",
                        default='hybpiper_stats')
    parser.add_argument('--hybpiper_dir',
                        default=None,
                        help='Specify directory containing HybPiper output sample folders. Default is the '
                             'current working directory.')
    parser.add_argument('--cpu',
                        type=int,
                        metavar='INTEGER',
                        help='Limit the number of CPUs. Default is to use all cores available minus '
                             'one.')

    parser.set_defaults(targetfile_dna=False, targetfile_aa=False)

    args = parser.parse_args()
    main(args)


def main(args):
    """
    Entry point for the hybpiper_main.py module.

    :param argparse.Namespace args:
    """

    logger.info(f'{"[INFO]:":10} HybPiper version {__version__} was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.info(f'{fill}\n')

    logger.info(f'{"[INFO]:":10} Recovering statistics for the HybPiper run(s)...')

    ####################################################################################################################
    # Set number of threads to use for steps that use multiprocessing:
    ####################################################################################################################
    if args.cpu:
        cpu = args.cpu
        logger.info(f'{"[INFO]:":10} Using {cpu} cpus/threads.')
    else:
        cpu = multiprocessing.cpu_count() - 1
        logger.info(f'{"[INFO]:":10} Number of cpus/threads not specified, using all available cpus minus 1 ({cpu}).')

    ####################################################################################################################
    # Set target file name and type:
    ####################################################################################################################
    targetfile = args.targetfile_dna if args.targetfile_dna else args.targetfile_aa
    targetfile_type = 'DNA' if args.targetfile_dna else 'protein'

    ####################################################################################################################
    # Check for presence of required input files:
    ####################################################################################################################
    logger.info(f'{"[INFO]:":10} The following file of sample names was provided: "{args.namelist}".')
    if not utils.file_exists_and_not_empty(args.namelist):
        sys.exit(f'{"[ERROR]:":10} File {args.namelist} is missing or empty, exiting!')

    logger.info(f'{"[INFO]:":10} The following target file was provided: "{targetfile}".')
    if not utils.file_exists_and_not_empty(targetfile):
        sys.exit(f'{"[ERROR]:":10} File {targetfile} is missing or empty, exiting!')

    ####################################################################################################################
    # Search within a user-supplied directory for the given sample directories, or the current directory if not:
    ####################################################################################################################
    if args.hybpiper_dir:
        sampledir_parent = os.path.abspath(args.hybpiper_dir)
    else:
        sampledir_parent = os.getcwd()

    if not os.path.isdir(sampledir_parent):
        sys.exit(f'{"[ERROR]:":10} Folder {sampledir_parent} not found, exiting!')

    ####################################################################################################################
    # Parse namelist and check for issues:
    ####################################################################################################################
    set_of_sample_names = utils.check_namelist(args.namelist,
                                               logger)

    ####################################################################################################################
    # Check if there is both an uncompressed folder AND a compressed file for any sample:
    ####################################################################################################################
    (samples_found,
     compressed_samples_set,
     uncompressed_samples_set) = utils.check_for_compressed_and_uncompressed_samples(set_of_sample_names,
                                                                                     sampledir_parent,
                                                                                     logger)

    ####################################################################################################################
    # Check for sample names present in the namelist.txt but not found in the directory provided:
    ####################################################################################################################
    list_of_sample_names = utils.check_for_missing_samples(args.namelist,
                                                           set_of_sample_names,
                                                           samples_found,
                                                           sampledir_parent,
                                                           logger)

    logger.info(f'{"[INFO]:":10} Total number of samples to process: {len(list_of_sample_names)}')

    ####################################################################################################################
    # Parse any compressed sample files using multiprocessing, to get contents and corresponding size dict:
    ####################################################################################################################
    compressed_sample_dict = dict()

    logger.info(f'{"[INFO]:":10} Getting compressed file contents for {len(compressed_samples_set)} compressed '
                f'samples...')

    with progressbar.ProgressBar(max_value=len(list_of_sample_names),
                                 min_poll_interval=30,
                                 widgets=widgets) as bar:

        with ProcessPoolExecutor(max_workers=cpu) as pool:
            manager = Manager()
            lock = manager.Lock()
            counter = manager.Value('i', 0)

            future_results = [pool.submit(utils.parse_compressed_sample,
                                          sample_name,
                                          sampledir_parent,
                                          lock,
                                          counter)

                              for sample_name in compressed_samples_set]

            for future in as_completed(future_results):
                try:
                    (sample_name,
                     sample_dict,
                     count) = future.result()

                    compressed_sample_dict[sample_name] = sample_dict

                    bar.update(count)

                except Exception as error:
                    print(f'Error raised: {error}')
                    tb = traceback.format_exc()
                    print(f'traceback is:\n{tb}')

            wait(future_results, return_when="ALL_COMPLETED")

    ####################################################################################################################
    # Get the names and lengths for each sequence in the target file:
    ####################################################################################################################
    gene_names = []
    reference_lengths = defaultdict(list)
    for ref_seq in SeqIO.parse(targetfile, "fasta"):
        ref_name = ref_seq.id.split("-")[-1]
        gene_names.append(ref_name)
        if targetfile_type.upper() == 'PROTEIN':
            reference_lengths[ref_name].append(len(ref_seq.seq) * 3)  # convert from amino-acids to nucleotides
        elif targetfile_type.upper() == 'DNA':
            reference_lengths[ref_name].append(len(ref_seq.seq))

    unique_gene_names = list(set(gene_names))
    avg_ref_lengths_dict = dict()
    for gene in sorted(unique_gene_names):
        avg_ref_lengths_dict[gene] = round(sum(reference_lengths[gene]) / len(reference_lengths[gene]))

    ####################################################################################################################
    # Iterate over each sample using multiprocessing:
    ####################################################################################################################

    logger.info(f'{"[INFO]:":10} Parsing data for all samples...')

    sample_stats_dict_collated = dict()
    warning_messages_dict = dict()

    with progressbar.ProgressBar(max_value=len(list_of_sample_names),
                                 min_poll_interval=30,
                                 widgets=widgets) as bar:

        with ProcessPoolExecutor(max_workers=cpu) as pool:
            manager = Manager()
            lock = manager.Lock()
            counter = manager.Value('i', 0)

            future_results = [pool.submit(parse_sample,
                                          sample_name,
                                          sampledir_parent,
                                          compressed_samples_set,
                                          compressed_sample_dict,
                                          unique_gene_names,
                                          args.sequence_type,
                                          lock,
                                          counter)

                              for sample_name in list_of_sample_names]

            for future in as_completed(future_results):
                try:
                    (sample_name,
                     sample_stats_dict,
                     count,
                     warning_messages_list) = future.result()

                    if sample_stats_dict:
                        sample_stats_dict_collated[sample_name] = sample_stats_dict

                    if len(warning_messages_list) != 0:
                        warning_messages_dict[sample_name] = warning_messages_list

                    bar.update(count)

                except Exception as error:
                    print(f'Error raised: {error}')
                    tb = traceback.format_exc()
                    print(f'traceback is:\n{tb}')

            wait(future_results, return_when="ALL_COMPLETED")

    ####################################################################################################################
    # Print any warning/info messages:
    ####################################################################################################################
    if len(warning_messages_dict) != 0:
        for sample, messages_list in sorted(warning_messages_dict.items(), key=lambda x: x[0]):
            logger.warning(f'{"[WARNING]:":10} For sample {sample}:')
            for message in messages_list:
                fill = utils.fill_forward_slash(f'{" "*10} {message}',
                                                width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                                break_on_forward_slash=True)
                logger.warning(fill)

    ####################################################################################################################
    # Write report file "seq_lengths.tsv":
    ####################################################################################################################
    seq_lengths_report_filename = f'{args.seq_lengths_filename}.tsv'
    lines_for_report = []

    # Get header and MeanLengths lines:
    sorted_ref_names_concat = '\t'.join([ref_name for ref_name in avg_ref_lengths_dict])  # already sorted above
    sorted_ref_mean_lengths_concat = '\t'.join([str(ref_mean_length) for ref_mean_length in avg_ref_lengths_dict.values()])
    header = f'Species\t{sorted_ref_names_concat}'
    mean_length_line = f'MeanLength\t{sorted_ref_mean_lengths_concat}'

    lines_for_report.append(header)
    lines_for_report.append(mean_length_line)

    # Get sample lines:
    for sample_name, stats_dict in sorted(sample_stats_dict_collated.items()):
        sample_seq_lengths = []
        for locus_name, locus_length in sorted(sample_stats_dict_collated[sample_name]['seq_lengths'].items(),
                                               key=lambda x: x[0]):
            ref_mean_length = avg_ref_lengths_dict[locus_name]

            if locus_length > 1.5 * ref_mean_length and args.sequence_type.upper() != 'SUPERCONTIG':
                logger.warning(f'{"[WARNING]:":10} Sequence length for {sample_name} is more than 50% '
                               f'longer than mean of references for {locus_name}!\n')
            sample_seq_lengths.append(str(locus_length))

        sample_seq_lengths_concat = '\t'.join(sample_seq_lengths)
        sample_line = f'{sample_name}\t{sample_seq_lengths_concat}'
        lines_for_report.append(sample_line)

    # Write the report:
    with open(seq_lengths_report_filename, 'w') as seq_lengths_handle:
        for item in lines_for_report:
            seq_lengths_handle.write(f'{item}\n')

    logger.info(f'{"[INFO]:":10} A sequence length table has been written to file: {args.seq_lengths_filename}.tsv')

    ####################################################################################################################
    # Write report file "hybpiper_stats.tsv":
    ####################################################################################################################
    seq_length_dict = seq_length_calc(seq_lengths_report_filename)

    lines_for_stats_report = []

    categories = [
        "Name",
        "NumReads",
        "ReadsMapped",
        "PctOnTarget",
        "GenesMapped",
        "GenesWithContigs",
        "GenesWithSeqs",
        "GenesAt25pct",
        "GenesAt50pct",
        "GenesAt75pct",
        "GenesAt150pct",
        "ParalogWarningsLong",
        "ParalogWarningsDepth",
        "GenesWithoutStitchedContigs",
        "GenesWithStitchedContigs",
        "GenesWithStitchedContigsSkipped",
        "GenesWithChimeraWarning",
        "TotalBasesRecovered"
    ]

    categories_for_printing = '\t'.join(categories)
    lines_for_stats_report.append(categories_for_printing)

    # Iterate over sample names and get recovery stats:
    stats_dict = {}

    for sample_name in sorted(sample_stats_dict_collated.keys()): # i.e. only samples that had stats recovered

        # Get the previously created stats dict for this sample:
        sample_stats_dict = sample_stats_dict_collated[sample_name]

        ################################################################################################################
        # Enrichment efficiency and recovery efficiency:
        ################################################################################################################
        stats_dict[sample_name] = [
            sample_stats_dict['total_input_reads'],
            sample_stats_dict['total_mapped_reads'],
            sample_stats_dict['pct_mapped'],
            sample_stats_dict['loci_with_mapped_reads'],
            sample_stats_dict['loci_with_contigs'],
            sample_stats_dict['loci_with_extracted_seqs'],
        ]

        ################################################################################################################
        # Loci with recovered lengths over percentage of mean refs length threshold stats:
        ################################################################################################################
        stats_dict[sample_name].extend(seq_length_dict[sample_name])

        ################################################################################################################
        # Paralog warning stats:
        ################################################################################################################
        stats_dict[sample_name].append(sample_stats_dict['paralog_warnings_long'])
        stats_dict[sample_name].append(sample_stats_dict['paralog_warnings_depth'])

        ################################################################################################################
        # Stitched contig stats:
        ################################################################################################################
        stats_dict[sample_name].append(sample_stats_dict['no_stitched_contig'])
        stats_dict[sample_name].append(sample_stats_dict['stitched_contig_produced'])
        stats_dict[sample_name].append(sample_stats_dict['stitched_contig_skipped'])

        ################################################################################################################
        # Stitched contig chimera warning stats:
        ################################################################################################################
        stats_dict[sample_name].append(sample_stats_dict['chimeric_contig_warnings'])

        ################################################################################################################
        # Total bases recovered from sample:
        ################################################################################################################
        stats_dict[sample_name].append(str(sample_stats_dict['sample_total_bases']))

    # Write to report file:
    for sample_name in stats_dict:
        stats_dict_for_writing = '\t'.join(stats_dict[sample_name])
        lines_for_stats_report.append(f'{sample_name}\t{stats_dict_for_writing}')

    with open(f'{args.stats_filename}.tsv', 'w') as hybpiper_stats_handle:
        for item in lines_for_stats_report:
            if len([item for item in item.split('\t')]) == 2:  # i.e. no bam file and no stats
                continue
            hybpiper_stats_handle.write(f'{item}\n')

    logger.info(f'{"[INFO]:":10} A statistics table has been written to file: {args.stats_filename}.tsv')
    logger.info(f'{"[INFO]:":10} Done!')


if __name__ == "__main__":
    standalone()
