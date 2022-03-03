#!/usr/bin/env python

"""
This script will retrieve paralog nucleotide (CDS) sequences for a specified gene in all samples located in
namelist.txt. It writes the unaligned sequences to folders, which can be user specified; defaults are 'paralogs_all'
and 'paralogs_no_chimeras'.

If a sample does not have paralogs for that gene, the sequence in the FNA directory is retrieved instead.
"""

import os
import sys
import argparse
from Bio import SeqIO
import logging


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
        if chimeric_genes_to_skip:
            logger.info(f'Putative chimeric gene sequences to skip from sample {sample_directory_name}:'
                        f' {chimeric_genes_to_skip}')
        else:
            logger.info(f'No putative chimeric gene sequences to skip from sample {sample_directory_name}')
    except FileNotFoundError:  # This file should be written in assemble.py even if it's empty
        logger.info(f'No chimeric stitched contig summary file found for gene sample {sample_directory_name}!')
        raise

    return chimeric_genes_to_skip


def retrieve_seqs(sample_base_directory_path, sample_directory_name, target_genes, fasta_dir_all=None, fasta_dir_no_chimeras=None):
    """
    Iterates over a list of gene name for a given sample, and for each gene produces two *.fasta files:  1) all
    paralog sequences (or *.FNA if no paralogs) from all samples; 2)  all paralog sequences (or non-chimeric *.FNA if
    no paralogs) from all samples.

    :param str sample_base_directory_path: path to the parent directory of the sample folder
    :param str sample_directory_name: name of the sample directory
    :param list target_genes: a list of target gene names
    :param str fasta_dir_all: folder name for fasta files with all paralogs
    :param str fasta_dir_no_chimeras: folder name for fasta files with all paralogs except putative chimeras
    :return list stats_for_report, genes_with_paralogs: lists of sequence counts and gene names for writing reports
    """

    # Make output directories:
    if not os.path.isdir(fasta_dir_all):
        os.mkdir(fasta_dir_all)
    if not os.path.isdir(fasta_dir_no_chimeras):
        os.mkdir(fasta_dir_no_chimeras)

    chimeric_genes_to_skip = get_chimeric_genes_for_sample(sample_directory_name)

    # # Recover a list of putative chimeric genes for the sample:
    # chimeric_genes_to_skip = []
    # try:
    #     with open(f'{sample_directory_name}/'
    #               f'{sample_directory_name}_genes_derived_from_putative_chimeric_stitched_contig.csv') as chimeric:
    #         lines = chimeric.readlines()
    #         for line in lines:
    #             chimeric_genes_to_skip.append(line.split(',')[1])
    #     if chimeric_genes_to_skip:
    #         logger.info(f'Putative chimeric gene sequences to skip from sample {sample_directory_name}:'
    #                     f' {chimeric_genes_to_skip}')
    #     else:
    #         logger.info(f'No putative chimeric gene sequences to skip from sample {sample_directory_name}')
    # except FileNotFoundError:  # This file should be written in assemble.py even if it's empty
    #     logger.info(f'No chimeric stitched contig summary file found for gene sample {sample_directory_name}!')
    #     raise

    # Normal recovery of all sequences; writes to folder fasta_dir_all:
    stats_for_report = [sample_directory_name]
    genes_with_paralogs = [sample_directory_name]
    seqs_to_write = None

    for gene in target_genes:
        has_paralogs = False

        # Get paths to paralogs file (might not be present if no paralogs) and *.FNA file:
        paralog_fasta_file_path = os.path.join(sample_base_directory_path,
                                               sample_directory_name,
                                               gene,
                                               sample_directory_name,
                                               'paralogs',
                                               f'{gene}_paralogs.fasta')

        fna_fasta_file_path = os.path.join(sample_base_directory_path,
                                           sample_directory_name,
                                           gene,
                                           sample_directory_name,
                                           'sequences',
                                           'FNA',
                                           f'{gene}.FNA')

        # Recover paralog sequences if present:
        if os.path.isfile(paralog_fasta_file_path):
            seqs_to_write = [x for x in SeqIO.parse(paralog_fasta_file_path, 'fasta')]
            num_seqs = len(seqs_to_write)
            has_paralogs = True

        # ...or recover nucleotide *.FNA sequence if present:
        elif os.path.isfile(fna_fasta_file_path):
            seqs_to_write = SeqIO.read(fna_fasta_file_path, 'fasta')
            num_seqs = 1

        if seqs_to_write:  # i.e. there were either paralogs or a *.FNA file for the gene
            SeqIO.write(seqs_to_write, f'{fasta_dir_all}/{sample_directory_name}_{gene}_paralogs_all.fasta', 'fasta')
        else:
            num_seqs = 0

        # Capture the number of sequences for this gene and add it to a list for report writing:
        stats_for_report.append(num_seqs)

        if has_paralogs:
            genes_with_paralogs.append(gene)

        # Skip any putative chimeric stitched contig sequences; writes to folder fasta_dir_no_chimeras:
        if gene in chimeric_genes_to_skip:
            logger.info(f'Skipping gene {gene} for sample {sample_directory_name}'
                        f' - putative chimeric stitched contig sequence!')
        else:
            if os.path.isfile(paralog_fasta_file_path):
                seqs_to_write = [x for x in SeqIO.parse(paralog_fasta_file_path, 'fasta')]

            elif os.path.isfile(fna_fasta_file_path):
                seqs_to_write = SeqIO.read(fna_fasta_file_path, 'fasta')

            if seqs_to_write:  # i.e. there were either paralogs or a *.FNA file for the gene
                SeqIO.write(seqs_to_write, f'{fasta_dir_no_chimeras}/'
                                           f'{sample_directory_name}_{gene}_paralogs_all.fasta', 'fasta')

    return stats_for_report, genes_with_paralogs
            

def standalone():
    """
    Used when this module is run as a stand-alone script. Parses command line arguments and runs function main().:
    """

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('namelist',
                        help='Text file containing list of HybPiper output directories, one per line.')
    parser.add_argument('targetfile',
                        help="FASTA file containing target sequences for each gene. Used to extract unique gene names "
                             "for paralog recovery")
    parser.add_argument('--fasta_dir_all',
                        help='Specify directory for output FASTA files (ALL)',
                        default='paralogs_all')
    parser.add_argument('--fasta_dir_no_chimeras',
                        help='Specify directory for output FASTA files (no putative chimeric sequences)',
                        default='paralogs_no_chimeras')
    parser.add_argument('--paralog_report_filename',
                        help='Specify the filename for the paralog *.tsv report table',
                        default='paralog_report')
    parser.add_argument('--genes_with_paralogs_filename',
                        help='Specify the filename for the *.txt list of genes with paralogs in at least one sample',
                        default='genes_with_paralogs')

    args = parser.parse_args()
    main(args)


def main(args):
    """
    Entry point for the assemble.py module.

    :param argparse.Namespace args:
    """

    # Get a list of genes to recover, parsed from the target file:
    target_genes = list(set([x.id.split('-')[-1] for x in SeqIO.parse(args.targetfile, 'fasta')]))

    # Get a list of sample names:
    namelist = [x.rstrip() for x in open(args.namelist)]

    stats_for_report_all_samples = []
    genes_with_paralogs_sample_lists = []
    for name in namelist:  # iterate over samples
        sample_base_directory_path, sample_directory_name = os.path.split(name)
        if not sample_directory_name:
            sample_base_directory_path, sample_folder_name = os.path.split(sample_base_directory_path)

        stats_for_report, genes_with_paralogs = retrieve_seqs(sample_base_directory_path,
                                                              sample_directory_name,
                                                              target_genes,
                                                              args.fasta_dir_all,
                                                              args.fasta_dir_no_chimeras)

        stats_for_report_all_samples.append(stats_for_report)
        genes_with_paralogs_sample_lists.append(genes_with_paralogs)

    # Get list of genes with paralogs in at least one sample, and write it to a report file:
    all_genes_with_paralogs_set = set()

    for gene_list in genes_with_paralogs_sample_lists:
        for item in gene_list[1:]:
            all_genes_with_paralogs_set.add(item)
    with open(f'{args.genes_with_paralogs_filename}.txt', 'w') as genes_with_paralogs_handle:
        for gene_name in all_genes_with_paralogs_set:
            genes_with_paralogs_handle.write(f'{gene_name}\n')

    # Write a *.tsv report, i.e. a matrix of sample names vs gene name with sequnce counts:
    with open(f'{args.paralog_report_filename}.tsv', 'w') as paralog_report_handle:
        genes = '\t'.join(target_genes)
        paralog_report_handle.write(f'\t{genes}\n')
        for stat_list in stats_for_report_all_samples:
            stats = '\t'.join([str(stat) for stat in stat_list])
            paralog_report_handle.write(f'{stats}\n')


if __name__ == "__main__":
    standalone()
