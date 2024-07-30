#!/usr/bin/env python

"""
This script will get the sequences generated from multiple runs of the 'hybpiper assemble' command.
Specify either a directory with all the HybPiper output directories or a file containing sample names of interest.
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


def get_chimeric_genes_for_sample(sampledir_parent,
                                  sample_name,
                                  compressed_sample_dict,
                                  compressed_sample_bool):
    """
    Returns a list of putative chimeric gene sequences for a given sample

    :param str sampledir_parent: directory name for the sample
    :param str sample_name: name of the sample
    :param dict compressed_sample_dict:
    :param bool compressed_sample_bool:
    :return list chimeric_genes_to_skip: a list of putative chimeric gene sequences for the sample
    """

    chimeric_genes_to_skip = []

    relative_path = f'{sample_name}/{sample_name}_genes_derived_from_putative_chimeric_stitched_contig.csv'
    full_path = (f'{sampledir_parent}/{sample_name}/'
                 f'{sample_name}_genes_derived_from_putative_chimeric_stitched_contig.csv')

    summary_file_not_found_bool = False

    if compressed_sample_bool:
        if relative_path in compressed_sample_dict[sample_name]:
            lines = utils.get_compressed_file_lines(sample_name,
                                                    sampledir_parent,
                                                    relative_path)
            for line in lines:
                chimeric_genes_to_skip.append(line.split(',')[1])
        else:
            summary_file_not_found_bool = True
    else:
        try:
            with open(full_path) as chimeric:
                lines = chimeric.readlines()
                for line in lines:
                    chimeric_genes_to_skip.append(line.split(',')[1])
        except FileNotFoundError:
            summary_file_not_found_bool = True

    if summary_file_not_found_bool:
        fill = textwrap.fill(f'{"[WARNING]:":10} No chimeric stitched contig summary file found for sample  '
                             f'{sample_name}. This usually occurs when no gene sequences were produced for '
                             f'this sample.', width=90, subsequent_indent=" " * 11)
        logger.warning(fill)

    return chimeric_genes_to_skip


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


def recover_sequences_from_all_samples(seq_dir,
                                       filename,
                                       target_genes,
                                       sample_names,
                                       hybpiper_dir=None,
                                       fasta_dir=None,
                                       skip_chimeric=False,
                                       stats_file=None,
                                       filter_by=False):
    """
    Recovers sequences (dna, amino acid, supercontig or intron) for all genes from all samples

    :param str seq_dir: directory to recover sequence from (FNA/FAA/intron)
    :param str filename: file name component used to reconstruct path to file (introns/supercontig/None)
    :param list target_genes: list of unique gene names in the target file
    :param str sample_names: text file with list of sample names
    :param None or str hybpiper_dir: if provided, a path to the directory containing HybPiper output
    :param None or str fasta_dir: directory name for output files, default is current directory
    :param bool skip_chimeric: if True, skip putative chimeric genes
    :param str stats_file: path to the stats file if provided
    :param list filter_by: a list of stats columns, 'greater'/'smaller' operators, and thresholds for filtering
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

    ####################################################################################################################
    # Search within a user-supplied directory for the given sample directories, or the current directory if not:
    ####################################################################################################################
    if hybpiper_dir:
        if os.path.isdir(hybpiper_dir):
            sampledir_parent = os.path.abspath(hybpiper_dir)
        else:
            sys.exit(f'Can not find a directory with name "{hybpiper_dir}, exiting..."')
    else:
        sampledir_parent = os.getcwd()

    ####################################################################################################################
    # Parse namelist and check for the presence of the corresponding sample directories or *.tar.gz files:
    ####################################################################################################################
    list_of_sample_names = []

    with open(sample_names, 'r') as namelist_handle:
        for line in namelist_handle.readlines():
            sample_name = line.rstrip()
            if sample_name:
                list_of_sample_names.append(sample_name)
                if re.search('/', sample_name):
                    sys.exit(f'{"[ERROR]:":10} A sample name must not contain '
                             f'forward slashes. The file {sample_names} contains: {sample_name}')

    samples_found = []
    compressed_sample_dict = dict()

    for sample_name in list_of_sample_names:
        compressed_sample = f'{sampledir_parent}/{sample_name}.tar.gz'
        uncompressed_sample = f'{sampledir_parent}/{sample_name}'

        if os.path.isfile(compressed_sample):
            compressed_sample_dict[sample_name] = utils.parse_compressed_sample(compressed_sample)
            samples_found.append(sample_name)

        if os.path.isdir(uncompressed_sample):
            if sample_name in samples_found:
                sys.exit(f'{"[ERROR]:":10} Both a compressed and an un-compressed sample folder have been found for '
                         f'sample {sample_name} in directory {sampledir_parent}. Please remove one!')
            else:
                samples_found.append(sample_name)

    samples_missing = sorted(list(set(list_of_sample_names) - set(samples_found)))

    if samples_missing:

        fill = utils.fill_forward_slash(f'{"[WARNING]:":10} File {sample_names} contains samples not found in '
                                        f'directory "{sampledir_parent}". The missing samples are:',
                                        width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                        break_on_forward_slash=True)

        logger.warning(f'{fill}\n')

        for name in samples_missing:
            logger.warning(f'{" " * 10} {name}')
        logger.warning('')

    list_of_sample_names = [x for x in list_of_sample_names if x not in samples_missing]

    ####################################################################################################################
    # Set output directory:
    ####################################################################################################################
    if fasta_dir:
        fasta_dir = os.path.abspath(fasta_dir)
        if not os.path.isdir(fasta_dir):
            os.mkdir(fasta_dir)
    else:
        fasta_dir = os.getcwd()

    if samples_to_recover:
        logger.info(f'{"[INFO]:":10} Retrieving {len(target_genes)} genes from {len(samples_to_recover)} samples')
    else:
        logger.info(f'{"[INFO]:":10} Retrieving {len(target_genes)} genes from {len(list_of_sample_names)} samples')

    ####################################################################################################################
    # Iterate over each gene:
    ####################################################################################################################
    samples_with_no_chimera_check_performed = set()
    for gene in target_genes:
        num_seqs = 0

        # Construct names for intron and supercontig output files:
        if seq_dir in ['intron', 'supercontig']:
            outfilename = f'{gene}_{filename}.fasta'
        else:
            outfilename = f'{gene}.{seq_dir}'

        with open(os.path.join(fasta_dir, outfilename), 'w') as outfile:

            ############################################################################################################
            # Iterate over each sample:
            ############################################################################################################
            for sample_name in list_of_sample_names:

                # Filter samples:
                if samples_to_recover and sample_name not in samples_to_recover:
                    continue

                # Check if the sample directory is a compressed tarball:
                compressed_sample_bool = True if sample_name in compressed_sample_dict else False

                ########################################################################################################
                # Determine whether a chimera check was performed for this sample during 'hybpiper assemble':
                ########################################################################################################
                chimera_check_performed_file = f'{sample_name}/{sample_name}_chimera_check_performed.txt'  # sample dir as root

                if compressed_sample_bool:
                    compressed_sample_bool_lines = utils.get_compressed_file_lines(sample_name,
                                                                                   sampledir_parent,
                                                                                   chimera_check_performed_file)
                    chimera_check_bool = compressed_sample_bool_lines[0]

                else:
                    with open(f'{sampledir_parent}/{chimera_check_performed_file}', 'r') as chimera_check_handle:
                        chimera_check_bool = chimera_check_handle.read()

                if chimera_check_bool == 'True':
                    chimera_check_performed_for_sample = True
                elif chimera_check_bool == 'False':
                    chimera_check_performed_for_sample = False
                else:
                    raise ValueError(f'chimera_check_bool is: {chimera_check_bool} for sample {sample_name}')

                ########################################################################################################
                # Recover a list of putative chimeric genes for the sample (if a chimera check was performed), and
                # skip gene if in the list:
                ########################################################################################################
                if not chimera_check_performed_for_sample:
                    samples_with_no_chimera_check_performed.add(sample_name)
                elif skip_chimeric:
                    chimeric_genes_to_skip = get_chimeric_genes_for_sample(sampledir_parent,
                                                                           sample_name,
                                                                           compressed_sample_dict,
                                                                           compressed_sample_bool)

                    if gene in chimeric_genes_to_skip:
                        logger.info(f'{"[INFO]:":10} Skipping putative chimeric stitched contig sequence for {gene}, '
                                    f'sample {sample_name}')
                        continue

                ########################################################################################################
                # Get path to the gene/intron/supercontig sequence with sample dir as root:
                ########################################################################################################
                if seq_dir == 'intron':
                    sample_path = os.path.join(sample_name, gene, sample_name, 'sequences', seq_dir, f'{gene}_{filename}.fasta')
                else:
                    sample_path = os.path.join(sample_name, gene, sample_name, 'sequences', seq_dir, f'{gene}.{seq_dir}')

                if compressed_sample_bool:
                    if sample_path in compressed_sample_dict[sample_name]:
                        sample_path_size = compressed_sample_dict[sample_name][sample_path]

                        if sample_path_size == 0:
                            logger.warning(f'{"[WARNING]:":10} File {sample_path} exists, but is empty!')
                        else:
                            seqrecord = utils.get_compressed_seqrecord(sample_name,
                                                                       sampledir_parent,
                                                                       sample_path)
                            SeqIO.write(seqrecord, outfile, 'fasta')
                            num_seqs += 1
                else:
                    full_sample_path = f'{sampledir_parent}/{sample_path}'
                    if os.path.isfile(full_sample_path):
                        if os.path.getsize(full_sample_path) == 0:
                            logger.warning(f'{"[WARNING]:":10} File {sample_path} exists, but is empty!')
                        else:
                            seqrecord = SeqIO.read(full_sample_path, 'fasta')
                            SeqIO.write(seqrecord, outfile, 'fasta')
                            num_seqs += 1

        logger.info(f'{"[INFO]:":10} Found {num_seqs} sequences for gene {gene}')

    # Warn user if --skip_chimeric_genes was provided but some samples didn't have a chimera check performed:
    if skip_chimeric and len(samples_with_no_chimera_check_performed) != 0:
        fill = textwrap.fill(f'{"[WARNING]:":10} Option "--skip_chimeric_genes" was provided but a chimera check '
                             f'was not performed during "hybpiper assemble" for the following samples:',
                             width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
        logger.warning(f'\n{fill}\n')

        for sample_name in sorted(list(samples_with_no_chimera_check_performed)):
            logger.warning(f'{" " * 10} {sample_name}')

        logger.warning(f'\n{" " * 10} No putative chimeric sequences were skipped for these samples!')


def recover_sequences_from_one_sample(seq_dir,
                                      filename,
                                      target_genes,
                                      single_sample_name,
                                      hybpiper_dir=None,
                                      fasta_dir=None,
                                      skip_chimeric=False):
    """
    Recovers sequences (dna, amino acid, supercontig or intron) for all genes from one sample

    :param str seq_dir: directory to recover sequence from
    :param str filename: file name component used to reconstruct path to file (None, intron or supercontig)
    :param list target_genes: list of unique gene names in the target file
    :param str single_sample_name: directory of a single sample
    :param None or str hybpiper_dir: if provided, a path to the directory containing HybPiper output
    :param None or str fasta_dir: directory name for output files, default is current directory
    :param bool skip_chimeric: if True, skip putative chimeric genes
    :return None:
    """

    ####################################################################################################################
    # Search within current directory or a user-supplied directory for the given sample directory:
    ####################################################################################################################
    if hybpiper_dir:
        if os.path.isdir(hybpiper_dir):
            sampledir_parent = os.path.abspath(hybpiper_dir)
        else:
            sys.exit(f'Can not find a directory with name "{hybpiper_dir}, exiting..."')
    else:
        sampledir_parent = os.getcwd()

    ####################################################################################################################
    # Check for the presence of the corresponding sample directory or *.tar.gz file:
    ####################################################################################################################
    sample_found = False
    compressed_sample_dict = dict()

    compressed_sample = f'{sampledir_parent}/{single_sample_name}.tar.gz'
    uncompressed_sample = f'{sampledir_parent}/{single_sample_name}'
    compressed_sample_bool = False

    if os.path.isfile(compressed_sample):
        compressed_sample_dict[single_sample_name] = utils.parse_compressed_sample(compressed_sample)
        sample_found = True
        compressed_sample_bool = True

    if os.path.isdir(uncompressed_sample):
        if sample_found:
            sys.exit(f'{"[ERROR]:":10} Both a compressed and an un-compressed sample folder have been found for '
                     f'sample {single_sample_name} in directory {sampledir_parent}. Please remove one!')
        else:
            sample_found = True

    if not sample_found:
        sys.exit(f'Can not find a directory or *.tar.gz file for sample "{single_sample_name}" within the directory '
                     f'"{sampledir_parent}", exiting...')

    ####################################################################################################################
    # Set output directory:
    ####################################################################################################################
    if fasta_dir:
        fasta_dir = os.path.abspath(fasta_dir)
        if not os.path.isdir(fasta_dir):
            os.mkdir(fasta_dir)
    else:
        fasta_dir = os.getcwd()

    logger.info(f'{"[INFO]:":10} Retrieving {len(target_genes)} genes from sample {single_sample_name}...')

    ####################################################################################################################
    # Construct names for intron and supercontig output files, and FNA/FAA files:
    ####################################################################################################################
    if seq_dir in ['intron', 'supercontig']:
        outfilename = f'{single_sample_name}_{filename}.fasta'
    else:
        outfilename = f'{single_sample_name}_{seq_dir}.fasta'

    ####################################################################################################################
    # Determine whether a chimera check was performed for this sample during 'hybpiper assemble':
    ####################################################################################################################
    chimera_check_performed_file = f'{single_sample_name}/{single_sample_name}_chimera_check_performed.txt'  # sample dir as root

    if compressed_sample_bool:
        compressed_sample_bool_lines = utils.get_compressed_file_lines(single_sample_name,
                                                                       sampledir_parent,
                                                                       chimera_check_performed_file)
        chimera_check_bool = compressed_sample_bool_lines[0]

    else:
        with open(f'{sampledir_parent}/{chimera_check_performed_file}', 'r') as chimera_check_handle:
            chimera_check_bool = chimera_check_handle.read()

    if chimera_check_bool == 'True':
        chimera_check_performed_for_sample = True
    elif chimera_check_bool == 'False':
        chimera_check_performed_for_sample = False
    else:
        raise ValueError(f'chimera_check_bool is: {chimera_check_bool} for sample {single_sample_name}')

    # Warn user if --skip_chimeric_genes was provided but some samples didn't have a chimera check performed:
    if skip_chimeric and not chimera_check_performed_for_sample:
        fill = textwrap.fill(f'{"[WARNING]:":10} Option "--skip_chimeric_genes" was provided but a chimera check '
                             f'was not performed during "hybpiper assemble" for sample {single_sample_name}. No '
                             f'putative chimeric sequences will be skipped for this sample!',
                             width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
        logger.warning(f'\n{fill}\n')

    ####################################################################################################################
    # Iterate over each gene:
    ####################################################################################################################
    sequences_to_write = []

    for gene in target_genes:
        num_seqs = 0

        # Recover a list of putative chimeric genes for the sample (if a chimera check was performed), and skip gene
        # if in list:
        if skip_chimeric and chimera_check_performed_for_sample:

            chimeric_genes_to_skip = get_chimeric_genes_for_sample(sampledir_parent,
                                                                   single_sample_name,
                                                                   compressed_sample_dict,
                                                                   compressed_sample_bool)

            if gene in chimeric_genes_to_skip:
                logger.info(f'{"[INFO]:":10} Skipping putative chimeric stitched contig sequence for {gene}, sample'
                            f' {single_sample_name}')
                continue

        ########################################################################################################
        # Get path to the gene/intron/supercontig sequence with sample dir as root:
        ########################################################################################################
        if seq_dir == 'intron':
            sample_path = os.path.join(single_sample_name, gene, single_sample_name, 'sequences', seq_dir,
                                       f'{gene}_{filename}.fasta')
        else:
            sample_path = os.path.join(single_sample_name, gene, single_sample_name, 'sequences', seq_dir,
                                       f'{gene}.{seq_dir}')

        if compressed_sample_bool:
            if sample_path in compressed_sample_dict[single_sample_name]:
                sample_path_size = compressed_sample_dict[single_sample_name][sample_path]

                if sample_path_size == 0:
                    logger.warning(f'{"[WARNING]:":10} File {sample_path} exists, but is empty!')
                else:
                    seqrecord = utils.get_compressed_seqrecord(single_sample_name,
                                                               sampledir_parent,
                                                               sample_path)
                    seqrecord.id = f'{seqrecord.id}-{gene}'
                    sequences_to_write.append(seqrecord)
                    num_seqs += 1
        else:
            full_sample_path = f'{sampledir_parent}/{sample_path}'
            if os.path.isfile(full_sample_path):
                if os.path.getsize(full_sample_path) == 0:
                    logger.warning(f'{"[WARNING]:":10} File {sample_path} exists, but is empty!')
                else:
                    seqrecord = SeqIO.read(full_sample_path, 'fasta')
                    seqrecord.id = f'{seqrecord.id}-{gene}'
                    sequences_to_write.append(seqrecord)
                    num_seqs += 1

        logger.info(f'{"[INFO]:":10} Found {num_seqs} sequences for gene {gene}')

    with open(os.path.join(fasta_dir, outfilename), 'w') as outfile:
        SeqIO.write(sequences_to_write, outfile, 'fasta')


def standalone():
    """
    Used when this module is run as a stand-alone script. Parses command line arguments and runs function main().:
    """

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    group_1 = parser.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna', dest='targetfile_dna', default=False,
                         help='FASTA file containing DNA target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa', dest='targetfile_aa', default=False,
                         help='FASTA file containing amino-acid target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    group_2 = parser.add_mutually_exclusive_group(required=True)
    group_2.add_argument('--sample_names',
                         help='Text file with names of HybPiper output directories, one per line.',
                         default=None)
    group_2.add_argument('--single_sample_name',
                         help='A single sample name to recover sequences for', default=None)
    parser.add_argument('sequence_type',
                        help='Type of sequence to extract',
                        choices=['dna', 'aa', 'intron', 'supercontig'])
    parser.add_argument('--hybpiper_dir', help='Specify directory containing HybPiper output',
                        default=None)
    parser.add_argument('--fasta_dir', help='Specify directory for output FASTA files',
                        default=None)
    parser.add_argument('--skip_chimeric_genes',
                        action='store_true',
                        dest='skip_chimeric',
                        help='Do not recover sequences for putative chimeric genes. This only has an effect for a '
                             'given sample if the option "--chimeric_stitched_contig_check" was provided to command '
                             '"hybpiper assemble".',
                        default=False)
    parser.add_argument('--stats_file', default=None,
                        help='Stats file produced by "hybpiper stats", required for selective filtering of retrieved '
                             'sequences')
    parser.add_argument('--filter_by', action='append', nargs=3,
                        help='Provide three space-separated arguments: 1) column of the stats_file to filter by, '
                             '2) "greater" or "smaller", 3) a threshold - either an integer (raw number of genes) or '
                             'float (percentage of genes in analysis).',
                        default=None)

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

    logger.info(f'{"[INFO]:":10} Recovering sequences for the HybPiper run(s)...')

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
    # Recover sequences from all samples:
    ####################################################################################################################
    if args.sample_names:
        recover_sequences_from_all_samples(seq_dir,
                                           filename,
                                           target_genes,
                                           args.sample_names,
                                           args.hybpiper_dir,
                                           args.fasta_dir,
                                           args.skip_chimeric,
                                           args.stats_file,
                                           args.filter_by)

    ####################################################################################################################
    # Recover sequences from a single sample:
    ####################################################################################################################
    elif args.single_sample_name:
        recover_sequences_from_one_sample(seq_dir,
                                          filename,
                                          target_genes,
                                          args.single_sample_name,
                                          args.hybpiper_dir,
                                          args.fasta_dir,
                                          args.skip_chimeric)


if __name__ == "__main__":
    standalone()
