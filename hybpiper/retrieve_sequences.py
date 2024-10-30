#!/usr/bin/env python

"""
This script will get the sequences generated from multiple runs of the 'hybpiper assemble' command.
Specify either a single sample name or a text file containing sample names of interest.
It retrieves all the gene names from the target file used in the run of the pipeline.

You must specify whether you want the protein (aa), nucleotide (dna) sequences.

You can also specify 'intron' to retrieve the intron sequences, or 'supercontig' to get intron and exon sequences.

Will output unaligned fasta files, one per gene.
"""

import os
import sys
import argparse
from Bio import SeqIO
import logging
import pandas
import textwrap
import progressbar
import multiprocessing
from multiprocessing import Manager
from concurrent.futures.process import ProcessPoolExecutor
from concurrent.futures import wait, as_completed
import traceback
import tarfile
from collections import defaultdict
import re

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


def get_samples_to_recover(filter_by, stats_df, target_genes):
    """
    Recovers a list of sample names that pass the filtering options requested. Returns the list

    :param Dataframe stats_df: pandas dataframe from the stats file
    :param list filter_by: a list of stats columns, 'greater'/'smaller' operators, and thresholds for filtering
    :param list target_genes: a list of unique target gene names in the targetfile
    :return list sample_to_retain: a list of sample names to retain after filtering
    """

    total_number_of_genes = len(target_genes)

    # Permanently changes the pandas settings
    # pandas.set_option('display.max_rows', None)
    # pandas.set_option('display.max_columns', None)
    # pandas.set_option('display.width', None)

    for filter_criterion in filter_by:
        column, operator, threshold = filter_criterion
        try:
            threshold = int(threshold)
            logger.info(f'{"[INFO]:":10} Threshold for {column} is: {operator} than {threshold}')
        except ValueError:
            fill = textwrap.fill(f'{"[INFO]:":10} Threshold {threshold} for {column} is a float: calculating as '
                                 'percentage of total number of genes in targetfile. Note that this threshold will be '
                                 'rounded down, if necessary.', width=90, subsequent_indent=" " * 11)
            logger.info(fill)
            threshold = round(float(threshold) * total_number_of_genes)
            logger.info(f'{"[INFO]:":10} Threshold for {column} is: {operator} than {threshold}')

        if operator == 'greater':
            stats_df = stats_df.loc[(stats_df[column] > threshold)]
        elif operator == 'smaller':
            stats_df = stats_df.loc[(stats_df[column] < threshold)]

    samples_to_retain = stats_df['Name'].tolist()

    return samples_to_retain


def parse_sample(sample_name,
                 sampledir_parent,
                 compressed_samples_set,
                 target_genes,
                 seq_dir,
                 filename,
                 lock,
                 counter):
    """

    :param str sample_name: name of the sample (no ".tar.gz")
    :param str sampledir_parent: path of the parent directory containing HybPiper output
    :param set compressed_samples_set: set of sample names that are compressed (i.e. *.tar.gz)
    :param list target_genes: list of unique gene names in the target file
    :param str seq_dir: directory to recover sequence from (FNA/FAA/intron)
    :param str filename: file name component used to reconstruct path to file (introns/supercontig/None)
    :param multiprocessing.managers.AcquirerProxy lock:
    :param multiprocessing.managers.ValueProxy counter:
    :return:
    """

    compressed_bool = True if sample_name in compressed_samples_set else False
    warning_messages_list = []
    sample_dict = dict()
    sample_dict['seqrecords'] = dict()

    ####################################################################################################################
    # Set paths of files to get/process for the sample:
    ####################################################################################################################
    locus_fasta_paths = []

    for gene_name in target_genes:
        if seq_dir == 'intron':
            locus_path = os.path.join(sample_name,
                                      gene_name,
                                      sample_name,
                                      'sequences',
                                      seq_dir,
                                      f'{gene_name}_{filename}.fasta')

            locus_fasta_paths.append(locus_path)
        else:
            locus_path = os.path.join(sample_name,
                                      gene_name,
                                      sample_name,
                                      'sequences',
                                      seq_dir,
                                      f'{gene_name}.{seq_dir}')

            locus_fasta_paths.append(locus_path)

    chimera_check_fn = f'{sample_name}/{sample_name}_chimera_check_performed.txt'
    chimeric_genes_list_fn = f'{sample_name}/{sample_name}_genes_derived_from_putative_chimeric_stitched_contig.csv'

    ####################################################################################################################
    # Set some default stats:
    ####################################################################################################################
    chimera_check = 'file_not_found'  # i.e. HybPiper version <2.2.0
    chimeric_genes_list = []

    ####################################################################################################################
    # Iterate over the files/folders for the sample and recover stats:
    ####################################################################################################################
    if compressed_bool:
        compressed_sample_file_path = f'{sampledir_parent}/{sample_name}.tar.gz'

        with tarfile.open(compressed_sample_file_path, 'r:gz') as tarfile_handle:

            for tarinfo in tarfile_handle.getmembers():
                if tarinfo.size != 0:

                    ####################################################################################################
                    # Get seqrecords:
                    ####################################################################################################
                    if tarinfo.name in locus_fasta_paths:
                        gene_name = os.path.basename(tarinfo.name).split('.')[0]
                        if re.search('_supercontig$', gene_name):
                            gene_name = gene_name.removesuffix('_supercontig')
                        elif re.search('_introns$', gene_name):
                            gene_name = gene_name.removesuffix('_introns')

                        sample_dict['seqrecords'][gene_name] = []

                        seqrecords = utils.get_compressed_seqrecords(tarfile_handle,
                                                                     tarinfo)
                        for seqrecord in seqrecords:
                            sample_dict['seqrecords'][gene_name].append(seqrecord)

                    ####################################################################################################
                    # Get chimera check bool (if present):
                    ####################################################################################################
                    if tarinfo.name == chimera_check_fn:
                        chimera_check_lines = utils.get_compressed_file_lines(tarfile_handle,
                                                                              tarinfo)
                        chimera_check = True if chimera_check_lines[0].strip() == 'True' else False

                    ####################################################################################################
                    # Get list of chimeric genes (if present):
                    ####################################################################################################
                    if tarinfo.name == chimeric_genes_list_fn:
                        chimeric_genes_list_lines = utils.get_compressed_file_lines(tarfile_handle,
                                                                                    tarinfo)
                        for line in chimeric_genes_list_lines:
                            chimeric_genes_list.append(line.split(',')[1])

    else:
        ################################################################################################################
        # Get seqrecords:
        ################################################################################################################
        for locus_path in locus_fasta_paths:
            locus_path_full = f'{sampledir_parent}/{locus_path}'
            gene_name = os.path.split(locus_path)[1].split('.')[0]
            if re.search('_supercontig$', gene_name):
                gene_name = gene_name.removesuffix('_supercontig')
            elif re.search('_introns$', gene_name):
                gene_name = gene_name.removesuffix('_introns')

            sample_dict['seqrecords'][gene_name] = []

            if utils.file_exists_and_not_empty(locus_path_full):
                seqrecords = SeqIO.parse(locus_path_full, 'fasta')

                for seqrecord in seqrecords:
                    sample_dict['seqrecords'][gene_name].append(seqrecord)

        ####################################################################################################
        # Get chimera check bool (if present):
        ####################################################################################################
        chimera_check_fn_full_path = f'{sampledir_parent}/{chimera_check_fn}'
        if utils.file_exists_and_not_empty(chimera_check_fn_full_path):

            with open(chimera_check_fn_full_path) as chimera_check_fn_handle:
                chimera_check_lines = list(chimera_check_fn_handle.readlines())
                chimera_check = True if chimera_check_lines[0].strip() == 'True' else False

        ####################################################################################################
        # Get list of chimeric genes (if present):
        ####################################################################################################
        chimeric_genes_list_fn_full_path = f'{sample_name}/{chimeric_genes_list_fn}'

        if utils.file_exists_and_not_empty(chimeric_genes_list_fn_full_path) == chimeric_genes_list_fn_full_path:

            with open(chimeric_genes_list_fn_full_path) as chimeric_genes_list_handle:
                for line in chimeric_genes_list_handle.readlines():
                    chimeric_genes_list.append(line.split(',')[1])

    ####################################################################################################################
    # Populate sample dict:
    ####################################################################################################################
    sample_dict['chimera_check'] = chimera_check
    sample_dict['chimeric_genes_list'] = chimeric_genes_list

    with lock:
        counter.value += 1
        return (sample_name,
                sample_dict,
                counter.value,
                warning_messages_list)


def recover_sequences_from_all_samples(list_of_sample_names,
                                       compressed_samples_set,
                                       seq_dir,
                                       filename,
                                       target_genes,
                                       sampledir_parent,
                                       fasta_dir,
                                       stats_file=None,
                                       filter_by=False,
                                       cpu=1,
                                       skip_chimeric_genes=False):
    """
    Recovers sequences (dna, amino acid, supercontig or intron) for all genes from all samples

    :param list list_of_sample_names: list of sample names to process (no missing samples)
    :param set compressed_samples_set: list of sample names that are *.tar.gz compressed files
    :param str seq_dir: directory to recover sequence from (FNA/FAA/intron)
    :param str filename: file name component used to reconstruct path to file (introns/supercontig/None)
    :param list target_genes: list of unique gene names in the target file
    :param None or str sampledir_parent: path of the parent directory containing HybPiper output
    :param None or str fasta_dir: directory name for output files, default is current directory
    :param str stats_file: path to the stats file if provided
    :param list filter_by: a list of stats columns, 'greater'/'smaller' operators, and thresholds for filtering
    :param int cpu: number of threads to use for multiprocessing
    :param bool skip_chimeric_genes: if True, skip putative chimeric genes
    :return None:
    """

    ####################################################################################################################
    # Read in the stats file, if present:
    ####################################################################################################################
    if stats_file:
        stats_df = pandas.read_csv(stats_file, delimiter='\t')
        samples_to_recover = get_samples_to_recover(filter_by, stats_df, target_genes)
        if not samples_to_recover:
            sys.exit(f'{"[ERROR]:":10} Your current filtering options will remove all samples! Please provide '
                     f'different filtering options!')
        logger.info(f'{"[INFO]:":10} The filtering options provided will recover sequences from'
                    f' {len(samples_to_recover)} sample(s). These are:')
        for sample_name in samples_to_recover:
            logger.info(f'{" " * 11}{sample_name}')
    else:
        samples_to_recover = False

    samples_to_parse = samples_to_recover if samples_to_recover else list_of_sample_names
    logger.info(f'{"[INFO]:":10} Retrieving {len(target_genes)} genes from {len(samples_to_parse)} samples')

    ####################################################################################################################
    # Iterate over each sample using multiprocessing:
    ####################################################################################################################

    logger.info(f'{"[INFO]:":10} Parsing data for all samples...')

    sample_dict_collated = dict()
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
                                          target_genes,
                                          seq_dir,
                                          filename,
                                          lock,
                                          counter)

                              for sample_name in list_of_sample_names]

            for future in as_completed(future_results):
                try:
                    (sample_name,
                     sample_dict,
                     count,
                     warning_messages_list) = future.result()

                    if sample_dict:
                        sample_dict_collated[sample_name] = sample_dict

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
    # Check for chimera_check value and provide warning message if not found for some samples:
    ####################################################################################################################
    samples_no_chimera_check_info = []

    for sample_name, sample_dict in sample_dict_collated.items():
        chimera_check = sample_dict['chimera_check']
        if chimera_check == 'file_not_found':
            samples_no_chimera_check_info.append(sample_name)

    if len(samples_no_chimera_check_info) != 0:
        chimera_file_fill = textwrap.fill(f'{"[WARNING]:":10} No file "<sample_name>_chimera_check_performed.txt" '
                                          f'found for one or more samples. This might be because you are running '
                                          f'"hybpiper retrieve_sequences" from Hybpiper version >=2.2.0 on a sample '
                                          f'that was assembled using HybPiper version <2.2.0. It will be assumed that '
                                          f'a chimera check was performed for this sample(s). Samples are:',
                                          width=90, subsequent_indent=" " * 11)

        logger.warning(chimera_file_fill)
        logger.info('')

        for sample_name in sorted(samples_no_chimera_check_info):
            logger.warning(f'{" " *10} {sample_name}')
            sample_dict_collated[sample_name]['chimera_check'] = True

    ####################################################################################################################
    # Warn user if --skip_chimeric_genes was provided but some samples didn't have a chimera check performed:
    ####################################################################################################################
    samples_with_no_chimera_check_performed = [sample_name for sample_name in sample_dict_collated if not
                                               sample_dict_collated[sample_name]['chimera_check']]

    if skip_chimeric_genes and len(samples_with_no_chimera_check_performed) != 0:
        fill = textwrap.fill(f'{"[WARNING]:":10} Option "--skip_chimeric_genes" was provided but a chimera check '
                             f'was not performed during "hybpiper assemble" for the following samples:',
                             width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
        logger.warning(f'\n{fill}\n')

        for sample_name in sorted(list(samples_with_no_chimera_check_performed)):
            logger.warning(f'{" " * 10} {sample_name}')

        logger.warning(f'\n{" " * 10} No putative chimeric sequences will be skipped for these samples!')

    ####################################################################################################################
    # Iterate over sample_dict_collated and collate loci seqrecords:
    ####################################################################################################################

    loci_to_write_dict = defaultdict(list)
    for sample_name, sample_dict in sample_dict_collated.items():
        seqrecord_dict = sample_dict['seqrecords']
        chimeric_genes_list = sample_dict['chimeric_genes_list']

        for locus_name, seqrecord_list in seqrecord_dict.items():
            if locus_name in chimeric_genes_list and skip_chimeric_genes:
                logger.info(f'{"[INFO]:":10} Skipping putative chimeric stitched contig sequence for {locus_name}, '
                            f'sample {sample_name}')
                continue
            else:
                for seqrecord in seqrecord_list:
                    loci_to_write_dict[locus_name].append(seqrecord)

    ####################################################################################################################
    # Write seqrecords to fasta file:
    ####################################################################################################################
    for locus_name, seqrecord_list in loci_to_write_dict.items():

        # Construct names for intron and supercontig output files:
        if seq_dir in ['intron', 'supercontig']:
            outfilename = f'{locus_name}.fasta'
        else:
            outfilename = f'{locus_name}.{seq_dir}'

        with open(os.path.join(fasta_dir, outfilename), 'w') as fasta_handle:
            SeqIO.write(seqrecord_list, fasta_handle, 'fasta')

        logger.info(f'{"[INFO]:":10} Found {len(seqrecord_list)} sequences for gene {locus_name}')

    logger.info(f'{"[INFO]:":10} Done!')


def recover_sequences_from_one_sample(seq_dir,
                                      filename,
                                      target_genes,
                                      single_sample_name,
                                      sampledir_parent,
                                      fasta_dir,
                                      skip_chimeric=False):
    """
    Recovers sequences (dna, amino acid, supercontig or intron) for all genes from one sample

    :param str seq_dir: directory to recover sequence from
    :param str filename: file name component used to reconstruct path to file (None, intron or supercontig)
    :param list target_genes: list of unique gene names in the target file
    :param str single_sample_name: directory of a single sample
    :param None or str sampledir_parent: path of the parent directory containing HybPiper output
    :param None or str fasta_dir: directory name for output files, default is current directory
    :param bool skip_chimeric: if True, skip putative chimeric genes
    :return None:
    """

    logger.info(f'{"[INFO]:":10} Retrieving {len(target_genes)} genes from sample {single_sample_name}...')

    ####################################################################################################################
    # Check for the presence of the corresponding sample directory or *.tar.gz file:
    ####################################################################################################################
    sample_found = False

    compressed_sample = f'{sampledir_parent}/{single_sample_name}.tar.gz'
    uncompressed_sample = f'{sampledir_parent}/{single_sample_name}'
    set_of_compressed_samples = set()

    if os.path.isfile(compressed_sample):
        sample_found = True
        set_of_compressed_samples.add(single_sample_name)

    if os.path.isdir(uncompressed_sample):
        if sample_found:
            sys.exit(f'{"[ERROR]:":10} Both a compressed and an un-compressed sample folder have been found for '
                     f'sample {single_sample_name} in directory {sampledir_parent}. Please remove one!')
        else:
            sample_found = True

    if not sample_found:
        sys.exit(f'Can not find a directory or "*.tar.gz" file for sample "{single_sample_name}" within the directory '
                 f'"{sampledir_parent}", exiting...')

    ####################################################################################################################
    # Process sample using multiprocessing (so parse_sample() can be used):
    ####################################################################################################################
    with ProcessPoolExecutor(max_workers=1) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)

        future_result = [pool.submit(parse_sample,
                                     single_sample_name,
                                     sampledir_parent,
                                     set_of_compressed_samples,
                                     target_genes,
                                     seq_dir,
                                     filename,
                                     lock,
                                     counter)]

        for future in as_completed(future_result):
            try:
                (sample_name,
                 sample_dict,
                 count,
                 warning_messages_list) = future.result()

            except Exception as error:
                print(f'Error raised: {error}')
                tb = traceback.format_exc()
                print(f'traceback is:\n{tb}')

        wait(future_result, return_when="ALL_COMPLETED")

    ####################################################################################################################
    # Print any warning/info messages:
    ####################################################################################################################
    if len(warning_messages_list) != 0:
        logger.warning(f'{"[WARNING]:":10} For sample {single_sample_name}:')
        for message in warning_messages_list:
            fill = utils.fill_forward_slash(f'{" " * 10} {message}',
                                            width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                            break_on_forward_slash=True)
            logger.warning(fill)

    ####################################################################################################################
    # Check for chimera_check value and exit with message if not found for this sample:
    ####################################################################################################################
    chimera_check = sample_dict['chimera_check']
    if chimera_check == 'file_not_found':
        chimera_file_fill = textwrap.fill(f'{"[WARNING]:":10} No file "<sample_name>_chimera_check_performed.txt" '
                                          f'found. This might be because you are running "hybpiper retrieve_sequences" '
                                          f'from Hybpiper version >=2.2.0 on a sample that was assembled using '
                                          f'HybPiper version <2.2.0. It will be assumed that a chimera check was '
                                          f'performed for this sample.',
                                          width=90, subsequent_indent=" " * 11)

        logger.warning(chimera_file_fill)

    ####################################################################################################################
    # Warn user if --skip_chimeric_genes was provided but some samples didn't have a chimera check performed:
    ####################################################################################################################
    if skip_chimeric and not chimera_check:
        fill = textwrap.fill(f'{"[WARNING]:":10} Option "--skip_chimeric_genes" was provided but a chimera check '
                             f'was not performed during "hybpiper assemble" for sample {single_sample_name}. No '
                             f'putative chimeric sequences will be skipped for this sample!',
                             width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
        logger.warning(f'\n{fill}\n')

    ####################################################################################################################
    # Parse sample_dict and collate loci seqrecords:
    ####################################################################################################################
    loci_to_write_dict = defaultdict(list)
    seqrecord_dict = sample_dict['seqrecords']
    chimeric_genes_list = sample_dict['chimeric_genes_list']

    for locus_name, seqrecord_list in seqrecord_dict.items():
        if locus_name in chimeric_genes_list and skip_chimeric:
            logger.info(f'{"[INFO]:":10} Skipping putative chimeric stitched contig sequence for {locus_name}, '
                        f'sample {sample_name}')
            continue
        else:
            for seqrecord in seqrecord_list:
                loci_to_write_dict[locus_name].append(seqrecord)

    ####################################################################################################################
    # Write seqrecords to fasta file:
    ####################################################################################################################
    sequences_to_write = []

    # Construct names for intron and supercontig output files:
    if seq_dir in ['intron', 'supercontig']:
        outfilename = f'{single_sample_name}_{filename}.fasta'
    else:
        outfilename = f'{single_sample_name}_{seq_dir}.fasta'

    for locus_name, seqrecord_list in loci_to_write_dict.items():

        for seqrecord in seqrecord_list:
            seqrecord.name = f'{seqrecord.name}-{locus_name}'
            seqrecord.id = f'{seqrecord.id}-{locus_name}'
            sequences_to_write.append(seqrecord)

        logger.info(f'{"[INFO]:":10} Found {len(seqrecord_list)} sequences for gene {locus_name}')

    with open(os.path.join(fasta_dir, outfilename), 'w') as fasta_handle:
        SeqIO.write(sequences_to_write, fasta_handle, 'fasta')

    logger.info(f'{"[INFO]:":10} Done!')


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
    group_2 = parser.add_mutually_exclusive_group(required=True)
    group_2.add_argument('--sample_names',
                         help='Text file with names of HybPiper output directories, one per line.',
                         default=None)
    group_2.add_argument('--single_sample_name',
                         help='A single sample name to recover sequences for',
                         default=None)
    parser.add_argument('sequence_type',
                        help='Type of sequence to extract',
                        choices=['dna', 'aa', 'intron', 'supercontig'])
    parser.add_argument('--hybpiper_dir',
                        help='Specify directory containing HybPiper output',
                        default=None)
    parser.add_argument('--fasta_dir',
                        help='Specify directory for output FASTA files',
                        default=None)
    parser.add_argument('--skip_chimeric_genes',
                        action='store_true',
                        dest='skip_chimeric',
                        help='Do not recover sequences for putative chimeric genes. This only has an effect for a '
                             'given sample if the option "--chimeric_stitched_contig_check" was provided to command '
                             '"hybpiper assemble".',
                        default=False)
    parser.add_argument('--stats_file',
                        default=None,
                        help='Stats file produced by "hybpiper stats", required for selective filtering of retrieved '
                             'sequences')
    parser.add_argument('--filter_by',
                        action='append',
                        nargs=3,
                        help='Provide three space-separated arguments: 1) column of the stats_file to filter by, '
                             '2) "greater" or "smaller", 3) a threshold - either an integer (raw number of genes) or '
                             'float (percentage of genes in analysis).',
                        default=None)
    parser.add_argument('--cpu',
                        type=int,
                        metavar='INTEGER',
                        help='Limit the number of CPUs. Default is to use all cores available minus '
                             'one.')

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

    logger.info(f'{"[INFO]:":10} Recovering {args.sequence_type} sequences from the HybPiper run(s)...')

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
    # Set target file name:
    ####################################################################################################################
    targetfile = args.targetfile_dna if args.targetfile_dna else args.targetfile_aa

    ####################################################################################################################
    # Check for presence of required input files:
    ####################################################################################################################
    if args.sample_names:
        logger.info(f'{"[INFO]:":10} The following file of sample names was provided: "{args.sample_names}".')
        if not utils.file_exists_and_not_empty(args.sample_names):
            sys.exit(f'{"[ERROR]:":10} File {args.sample_names} is missing or empty, exiting!')

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
    # Set output directory:
    ####################################################################################################################
    if args.fasta_dir:
        fasta_dir = os.path.abspath(args.fasta_dir)
        if not os.path.isdir(fasta_dir):
            os.mkdir(fasta_dir)
    else:
        fasta_dir = os.getcwd()

    ####################################################################################################################
    # Check some args:
    ####################################################################################################################
    columns = ['GenesMapped',
               'GenesWithContigs',
               'GenesWithSeqs',
               'GenesAt25pct',
               'GenesAt50pct',
               'GenesAt75pct',
               'GenesAt150pct',
               'ParalogWarningsLong',
               'ParalogWarningsDepth',
               'GenesWithoutSupercontigs',
               'GenesWithSupercontigs',
               'GenesWithSupercontigSkipped',
               'GenesWithChimeraWarning',
               'TotalBasesRecovered']

    operators = ['greater', 'smaller']

    if args.stats_file and not args.filter_by:
        sys.exit(f'{"[ERROR]:":10} A stats file has been provided but no filtering options have been specified via '
                 f'the parameter "--filter_by". Either provide filtering criteria or omit the "--stats_file" '
                 f'parameter!')

    if args.filter_by and not args.stats_file:
        sys.exit(f'{"[ERROR]:":10} Filtering options has been provided but no stats file has been specified via '
                 f'the parameter "--stats_file". Either provide the stats file or omit the "--filter_by" parameter(s)!')

    if args.filter_by:
        columns_to_filter = [item[0] for item in args.filter_by]
        operators_to_filter = [item[1] for item in args.filter_by]
        thresholds_to_filter = [item[2] for item in args.filter_by]

        if not all(column in columns for column in columns_to_filter):
            sys.exit(f'{"[ERROR]:":10} Only columns from the following list are allowed: {columns}')
        if not all(operator in operators for operator in operators_to_filter):
            sys.exit(f'{"[ERROR]:":10} Only operators from the following list are allowed: {operators}')

        for threshold in thresholds_to_filter:
            try:
                threshold_is_float = float(threshold)
            except ValueError:
                sys.exit(f'{"[ERROR]:":10} Please provide only integers or floats as threshold values. You have '
                         f'provided: {threshold}')

    ####################################################################################################################
    # Set sequence directory name and file names:
    ####################################################################################################################
    seq_dir = None
    filename = None

    if args.sequence_type == 'dna':
        seq_dir = "FNA"
    elif args.sequence_type == 'aa':
        seq_dir = "FAA"
    elif args.sequence_type == 'intron':
        seq_dir = 'intron'
        filename = 'introns'
    elif args.sequence_type == 'supercontig':
        seq_dir = 'intron'
        filename = 'supercontig'

    assert seq_dir

    ####################################################################################################################
    # Get gene names parsed from a target file.
    ####################################################################################################################
    target_genes = list(set([x.id.split('-')[-1] for x in SeqIO.parse(targetfile, 'fasta')]))

    ####################################################################################################################
    # Parse namelist and check for issues:
    ####################################################################################################################
    if args.sample_names:  # i.e. it's not for a single sample
        set_of_sample_names = utils.check_namelist(args.sample_names,
                                                   logger)

        ################################################################################################################
        # Check if there is both an uncompressed folder AND a compressed file for any sample:
        ################################################################################################################
        (samples_found,
         compressed_samples_set,
         uncompressed_samples_set) = utils.check_for_compressed_and_uncompressed_samples(set_of_sample_names,
                                                                                         sampledir_parent,
                                                                                         logger)

        ################################################################################################################
        # Check for sample names present in the namelist.txt but not found in the directory provided:
        ################################################################################################################
        list_of_sample_names = utils.check_for_missing_samples(args.sample_names,
                                                               set_of_sample_names,
                                                               samples_found,
                                                               sampledir_parent,
                                                               logger)

        ################################################################################################################
        # Recover sequences from all samples:
        ################################################################################################################
        recover_sequences_from_all_samples(list_of_sample_names,
                                           compressed_samples_set,
                                           seq_dir,
                                           filename,
                                           target_genes,
                                           sampledir_parent,
                                           fasta_dir,
                                           args.stats_file,
                                           args.filter_by,
                                           cpu,
                                           args.skip_chimeric,)

    ###################################################################################################################
    # Recover sequences from a single sample:
    ###################################################################################################################
    elif args.single_sample_name:
        recover_sequences_from_one_sample(seq_dir,
                                          filename,
                                          target_genes,
                                          args.single_sample_name,
                                          sampledir_parent,
                                          fasta_dir,
                                          args.skip_chimeric)


if __name__ == "__main__":
    standalone()
