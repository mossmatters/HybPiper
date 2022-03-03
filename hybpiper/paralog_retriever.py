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
    :return xxx num_seqs: xxx
    """

    chimeric_genes_to_skip = []
    seqs_to_write = None

    # Make output directories:
    if not os.path.isdir(fasta_dir_all):
        os.mkdir(fasta_dir_all)
    if not os.path.isdir(fasta_dir_no_chimeras):
        os.mkdir(fasta_dir_no_chimeras)

    # Recover a list of putative chimeric genes for the sample:
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

    # Normal recovery of all sequences; writes to folder fasta_dir_all:
    stats_for_report = [sample_directory_name]

    for gene in target_genes:

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

    return stats_for_report
            

def standalone():
    """
    Used when this module is run as a stand-alone script. Parses command line arguments and runs function main().:
    """

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('namelist', help='Text file containing list of HybPiper output directories, one per line.')
    parser.add_argument('targetfile', help="FASTA file containing target sequences for each gene. Used to extract "
                                           "unique gene names for paralog recovery")
    parser.add_argument('--fasta_dir_all', help='Specify directory for output FASTA files (ALL)',
                        default='paralogs_all')
    parser.add_argument('--fasta_dir_no_chimeras', help='Specify directory for output FASTA files (no putative '
                                                        'chimeric sequences)', default='paralogs_no_chimeras')
    parser.add_argument('--paralog_report_filename', help='Specify the filename for the paralog *.tsv report table',
                        default='paralog_report')

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
    for name in namelist:  # iterate over samples
        sample_base_directory_path, sample_directory_name = os.path.split(name)
        if not sample_directory_name:
            sample_base_directory_path, sample_folder_name = os.path.split(sample_base_directory_path)

        stats_for_report_all_samples.append(retrieve_seqs(sample_base_directory_path,
                                            sample_directory_name,
                                            target_genes,
                                            args.fasta_dir_all,
                                            args.fasta_dir_no_chimeras))

    # Write a *.tsv report:
    with open(f'{args.paralog_report_filename}.tsv', 'w') as paralog_report:
        genes = '\t'.join(target_genes)
        paralog_report.write(f'\t{genes}\n')
        for stat_list in stats_for_report_all_samples:
            stats = '\t'.join([str(stat) for stat in stat_list])
            paralog_report.write(f'{stats}\n')


if __name__ == "__main__":
    standalone()
