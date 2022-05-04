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


# Create a custom logger

# Log to Terminal (stderr):
console_handler = logging.StreamHandler(sys.stderr)
console_handler.setLevel(logging.INFO)

# Setup logger:
logger = logging.getLogger(f'hybpiper.{__name__}')

# Add handlers to the logger
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)  # Default level is 'WARNING'


def get_chimeric_genes_for_sample(sample_directory_name):
    """
    Returns a list of putative chimeric gene sequences for a given sample

    :param str sample_directory_name: directory name for the sample
    :return list chimeric_genes_to_skip: a list of putative chimeric gene sequences for the sample
    """

    chimeric_genes_to_skip = []
    try:
        with open(f'{sample_directory_name}/'
                  f'{sample_directory_name}_genes_derived_from_putative_chimeric_stitched_contig.csv') as chimeric:
            lines = chimeric.readlines()
            for line in lines:
                chimeric_genes_to_skip.append(line.split(',')[1])
    except FileNotFoundError:  # This file should be written in assemble.py even if it's empty
        logger.info(f'No chimeric stitched contig summary file found for sample {sample_directory_name}!')
        raise

    return chimeric_genes_to_skip


def get_samples_to_recover(filter_by, stats_df, target_genes):
    """
    Recovers a list of sample names that pass the filtering options requested. Returns the list

    :param str stats_df: pandas dataframe from the stats file
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


def recover_sequences_from_all_samples(seq_dir, filename, target_genes, sample_names, hybpiper_dir=None,
                                       fasta_dir=None, skip_chimeric=False, stats_file=None, filter_by=False):
    """
    Recovers sequences (dna, amino acid, supercontig or intron) for all genes from all samples

    :param str seq_dir: directory to recover sequence from
    :param str filename: file name component used to reconstruct path to file
    :param list target_genes: list of unique gene names in the target file
    :param str sample_names: directory of samples, or text file with list of sample names
    :param None or str hybpiper_dir: if provided, a path to the directory containing HybPiper output
    :param None or str fasta_dir: directory name for output files, default is current directory
    :param bool skip_chimeric: if True, skip putative chimeric genes
    :param str stats_file: path to the stats file if provided
    :param list filter_by: a list of stats columns, 'greater'/'smaller' operators, and thresholds for filtering
    :return None:
    """

    # Read in the stats file, if present:
    if stats_file:
        stats_df = pandas.read_csv(stats_file, delimiter='\t')
        samples_to_recover = get_samples_to_recover(filter_by, stats_df, target_genes)
        if not samples_to_recover:
            sys.exit(f'{"[ERROR]:":10} Your current filtering options will remove all samples! Please provide '
                     f'different filtering options!')
        logger.info(f'{"[INFO]:":10} The filtering options provided will recover sequences from'
                    f' {len(samples_to_recover)} sample(s). These are:')
        for sample in samples_to_recover:
            logger.info(f'{" " * 11}{sample}')
    else:
        samples_to_recover = False

    # Recover sample names from directory names or names in text file provided:
    if os.path.isdir(sample_names):
        sampledir = sample_names
        sample_names = [x for x in os.listdir(sampledir) if os.path.isdir(os.path.join(sampledir, x)) and not
                        x.startswith('.')]
    else:
        sample_names = [x.rstrip() for x in open(sample_names)]
        if hybpiper_dir:
            sampledir = hybpiper_dir
        else:
            sampledir = '.'

    # Set output directory:
    if fasta_dir:
        fasta_dir = fasta_dir
        if not os.path.isdir(fasta_dir):
            os.mkdir(fasta_dir)
    else:
        fasta_dir = '.'

    if samples_to_recover:
        logger.info(f'{"[INFO]:":10} Retrieving {len(target_genes)} genes from {len(samples_to_recover)} samples')
    else:
        logger.info(f'{"[INFO]:":10} Retrieving {len(target_genes)} genes from {len(sample_names)} samples')
    for gene in target_genes:
        numSeqs = 0

        # Construct names for intron and supercontig output files:
        if seq_dir in ['intron', 'supercontig']:
            outfilename = f'{gene}_{filename}.fasta'
        else:
            outfilename = f'{gene}.{seq_dir}'

        with open(os.path.join(fasta_dir, outfilename), 'w') as outfile:
            for sample in sample_names:

                # Filter samples:
                if samples_to_recover and sample not in samples_to_recover:
                    # print(f'Sample {sample} did not pass filtering criteria, skipping recovery...')
                    continue

                # Recover a list of putative chimeric genes for the sample, and skip gene if in list:
                if skip_chimeric:
                    chimeric_genes_to_skip = get_chimeric_genes_for_sample(sample)
                    # print(f'chimeric_genes_to_skip is: {chimeric_genes_to_skip}')
                    if gene in chimeric_genes_to_skip:
                        logger.info(f'Skipping putative chimeric stitched contig sequence for {gene}, sample {sample}')
                        continue

                # Get path to the gene/intron/supercontig sequence:
                if seq_dir == 'intron':
                    sample_path = os.path.join(sampledir, sample, gene, sample, 'sequences', seq_dir,
                                               f'{gene}_{filename}.fasta')
                else:
                    sample_path = os.path.join(sampledir, sample, gene, sample, 'sequences', seq_dir,
                                               f'{gene}.{seq_dir}')
                try:
                    seq = next(SeqIO.parse(sample_path, 'fasta'))
                    # print(seq)
                    SeqIO.write(seq, outfile, 'fasta')
                    numSeqs += 1
                # except FileNotFoundError or StopIteration:  # BioPython 1.80 returns StopIteration error?
                except FileNotFoundError:
                    pass
        logger.info(f'{"[INFO]:":10} Found {numSeqs} sequences for gene {gene}')


def recover_sequences_from_one_sample(seq_dir,
                                      filename,
                                      target_genes,
                                      single_sample_name,
                                      fasta_dir=None,
                                      skip_chimeric=False):
    """
    Recovers sequences (dna, amino acid, supercontig or intron) for all genes from one sample

    :param str seq_dir: directory to recover sequence from
    :param str filename: file name component used to reconstruct path to file (None, intron or supercontig)
    :param list target_genes: list of unique gene names in the target file
    :param str single_sample_name: directory of a single sample
    :param None or str fasta_dir: directory name for output files, default is current directory
    :param bool skip_chimeric: if True, skip putative chimeric genes
    :return None:
    """

    if not os.path.isdir(single_sample_name):
        sys.exit(f'Can not find a directory for sample {single_sample_name}, exiting...')

    # Create a user-supplied directory if provided, or write to the current directory if not:
    if fasta_dir:
        fasta_dir = fasta_dir
        if not os.path.isdir(fasta_dir):
            os.mkdir(fasta_dir)
    else:
        fasta_dir = '.'
    logger.info(f'{"[INFO]:":10} Retrieving {len(target_genes)} genes from sample {single_sample_name}...')
    sequences_to_write = []

    # Construct names for intron and supercontig output files:
    if seq_dir in ['intron', 'supercontig']:
        outfilename = f'{filename}.fasta'
    else:
        outfilename = f'{seq_dir}.fasta'
    for gene in target_genes:
        numSeqs = 0
        # Recover a list of putative chimeric genes for the sample, and skip gene if in list:
        if skip_chimeric:
            chimeric_genes_to_skip = get_chimeric_genes_for_sample(single_sample_name)
            # print(f'chimeric_genes_to_skip is: {chimeric_genes_to_skip}')
            if gene in chimeric_genes_to_skip:
                logger.info(f'{"[INFO]:":10} Skipping putative chimeric stitched contig sequence for {gene}, sample'
                            f' {single_sample_name}')
                continue

        # Get path to the gene/intron/supercontig sequence:
        if seq_dir == 'intron':
            sample_path = os.path.join(single_sample_name, gene, single_sample_name, 'sequences', seq_dir,
                                       f'{gene}_{filename}.fasta')
        else:
            sample_path = os.path.join(single_sample_name, gene, single_sample_name, 'sequences', seq_dir,
                                       f'{gene}.{seq_dir}')
        try:
            seq = next(SeqIO.parse(sample_path, 'fasta'))
            seq.id = f'{seq.id}-{gene}'
            sequences_to_write.append(seq)
            numSeqs += 1
        # except FileNotFoundError or StopIteration:  # BioPython 1.80 returns StopIteration error?
        except FileNotFoundError:
            pass
        logger.info(f'Found {numSeqs} sequences for gene {gene}.')

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
                         help='Directory containing Hybpiper output OR a file containing HybPiper output names, '
                              'one per line',
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
                        help='Do not recover sequences for putative chimeric genes',
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
    Entry point for the assemble.py module.

    :param argparse.Namespace args:
    """

    # Check some args:
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
               'GenesWithChimeraWarning']

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
            sys.exit(f'Only columns from the following list are allowed: {columns}')
        if not all(operator in operators for operator in operators_to_filter):
            sys.exit(f'Only operators from the following list are allowed: {operators}')
        for threshold in thresholds_to_filter:
            try:
                threshold_is_float = float(threshold)
            except ValueError:
                sys.exit(f'Please provide only integers or floats as threshold values. You have provided: {threshold}')

    # Set target file name:
    if args.targetfile_dna:
        targetfile = args.targetfile_dna
    elif args.targetfile_aa:
        targetfile = args.targetfile_aa

    # Set sequence directory name and file names:
    if args.sequence_type == 'dna':
        seq_dir = "FNA"
        filename = None
    elif args.sequence_type == 'aa':
        seq_dir = "FAA"
        filename = None
    elif args.sequence_type == 'intron':
        seq_dir = 'intron'
        filename = 'introns'
    elif args.sequence_type == 'supercontig':
        seq_dir = 'intron'
        filename = 'supercontig'

    # Use gene names parsed from a target file.
    target_genes = list(set([x.id.split('-')[-1] for x in SeqIO.parse(targetfile, 'fasta')]))

    # Recover sequences from all samples:
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
    elif args.single_sample_name:
        recover_sequences_from_one_sample(seq_dir,
                                          filename,
                                          target_genes,
                                          args.single_sample_name,
                                          args.fasta_dir,
                                          args.skip_chimeric)


if __name__ == "__main__":
    standalone()
