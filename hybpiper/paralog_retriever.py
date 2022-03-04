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
from collections import defaultdict

try:
    import pandas as pd
except ImportError:
    sys.exit(f"Required Python package 'pandas' not found. Is it installed for the Python used to run this script?")

try:
    import seaborn as sns
except ImportError:
    sys.exit(f"Required Python package 'seaborn' not found. Is it installed for the Python used to run this script?")

try:
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit(f"Required Python package 'matplotlib' not found. Is it installed for the Python used to run this script?")

# Create a custom logger

# Log to Terminal (stderr):
console_handler = logging.StreamHandler(sys.stderr)
console_handler.setLevel(logging.INFO)

# Setup logger:
logger = logging.getLogger(f'hybpiper.{__name__}')

# Add handlers to the logger
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)  # Default level is 'WARNING'


def get_figure_dimensions(df, figure_length, figure_height, sample_text_size, gene_text_size):
    """
    Takes a dataframe and returns figure length (inches), figure height (inches), sample_text_size and gene_id_text_size
    values based on the number of samples and genes in the paralog_report.tsv input file provided.

    :param pandas.core.frame.DataFrame df: pandas dataframe of paralog_report.tsv after filtering and pivot
    :param NoneType or int figure_length: if provided, dimension (in inches) for the figure length
    :param NoneType or int figure_height: if provided, dimension (in inches) for the figure height
    :param NoneType or int sample_text_size: if provided, dimension (in inches) for the figure sample text size
    :param NoneType or int gene_text_size: if provided, dimension (in inches) for the figure gene_id text size
    :return float fig_length, figure_height, sample_text_size, gene_id_text_size:
    """

    num_samples = len(df.index)
    num_genes = len(df.columns)

    logger.info(f'Number of samples in input paralog file is: {num_samples}')
    logger.info(f'Number of genes in input paralog file is: {num_genes}')

    # Set default text label size (in points) unless specified at the command line:
    sample_text_size = sample_text_size if sample_text_size else 10
    gene_id_text_size = gene_text_size if gene_text_size else 10

    # Set figure height dimensions for a given number of samples:
    if num_samples <= 10:
        figure_height = figure_height if figure_height else 4
    elif 10 < num_samples <= 20:
        figure_height = figure_height if figure_height else 8
    elif 20 < num_samples <= 50:
        figure_height = figure_height if figure_height else 16
    elif 50 < num_samples <= 100:
        figure_height = figure_height if figure_height else 16
    elif 100 < num_samples <= 200:
        figure_height = figure_height if figure_height else 140
    elif 200 < num_samples <= 400:
        figure_height = figure_height if figure_height else 180
    elif num_samples > 400:
        figure_height = figure_height if figure_height else 240

    # Set figure length dimensions for a given number of genes (i.e. number of unique genes in target file):
    if num_genes <= 10:
        fig_length = figure_length if figure_length else 4
    elif 10 < num_genes <= 20:
        fig_length = figure_length if figure_length else 6
    elif 20 < num_genes <= 50:
        fig_length = figure_length if figure_length else 8
    elif 50 < num_genes <= 100:
        fig_length = figure_length if figure_length else 8
    elif 100 < num_genes <= 200:
        fig_length = figure_length if figure_length else 70
    elif 200 < num_genes <= 400:
        fig_length = figure_length if figure_length else 90
    elif num_genes > 400:
        fig_length = figure_length if figure_length else 120

    logger.info(f'fig_length: {fig_length} inches, figure_height: {figure_height} inches, sample_text_size:'
                f' {sample_text_size} points, gene_id_text_size: {gene_id_text_size} points')

    return fig_length, figure_height, sample_text_size, gene_id_text_size


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


def retrieve_seqs(sample_base_directory_path, sample_directory_name, target_genes, fasta_dir_all=None,
                  fasta_dir_no_chimeras=None):
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


def write_paralogs_above_threshold_report(genes_with_paralogs_sample_lists, paralogs_list_threshold_percentage,
                                          namelist, target_genes, paralogs_above_threshold_report_filename):
    """
    Writes a *.txt report list gene names and sample names for genes and samples that have paralogs in >=

    :param list genes_with_paralogs_sample_lists: list of lists containing ['sample_name', gene_with_paralog1, 'etc']
    :param float paralogs_list_threshold_percentage: percentage threshold used for report calculations
    :param list namelist: list of sample names
    :param list target_genes: list of target gen names
    :param str paralogs_above_threshold_report_filename: file name for the report file
    :return:
    """

    # Calculate number of genes with warnings in >= paralogs_list_threshold_percentage samples:
    genes_with_paralogs_in_greater_than_threshold_samples = []
    gene_paralog_count_dict = defaultdict(int)

    for sample_gene_paralog_warning_list in genes_with_paralogs_sample_lists:
        for gene_name in sample_gene_paralog_warning_list[1:]:
            # print(gene_name)
            gene_paralog_count_dict[gene_name] += 1
    print(f'gene_paralog_count_dict is: {gene_paralog_count_dict}')

    for gene, paralog_count in gene_paralog_count_dict.items():
        if paralog_count / len(namelist) >= paralogs_list_threshold_percentage:  # Check percentage of samples
            # print(gene)
            # print(paralog_count/len(namelist))
            genes_with_paralogs_in_greater_than_threshold_samples.append(gene)

    logger.info(f'There are {len(genes_with_paralogs_in_greater_than_threshold_samples)} genes with paralogs in >= '
                f'{paralogs_list_threshold_percentage}% of samples. The gene names have been written to the '
                f'file "{paralogs_above_threshold_report_filename}.txt"')

    # Calculate number of samples with warnings in >= paralogs_list_threshold_percentage genes:
    samples_with_paralogs_in_greater_than_threshold_genes = []
    sample_paralog_count_dict = defaultdict(int)

    for sample_gene_paralog_warning_list in genes_with_paralogs_sample_lists:
        # print(sample_gene_paralog_warning_list)
        sample_name = sample_gene_paralog_warning_list[0]
        for gene_name in sample_gene_paralog_warning_list[1:]:
            # print(gene_name)
            sample_paralog_count_dict[sample_name] += 1
    print(f'sample_paralog_count_dict is: {sample_paralog_count_dict}')

    for sample, paralog_count in sample_paralog_count_dict.items():
        if paralog_count / len(target_genes) >= paralogs_list_threshold_percentage:
            # print(sample)
            # print(paralog_count/len(target_genes))
            samples_with_paralogs_in_greater_than_threshold_genes.append(sample)

    logger.info(f'There are {len(samples_with_paralogs_in_greater_than_threshold_genes)} samples with paralogs in >= '
                f'{paralogs_list_threshold_percentage}% of genes. The sample names have been written to the '
                f'file "{paralogs_above_threshold_report_filename}.txt"')

    # Write the text report file:
    with open(f'{paralogs_above_threshold_report_filename}.txt', 'w') as genes_with_paralogs_handle:
        genes_with_paralogs_handle.write(f'There are {len(samples_with_paralogs_in_greater_than_threshold_genes)} '
                                         f'samples with paralogs in >= {paralogs_list_threshold_percentage}% of '
                                         f'genes. The sample names are:\n')
        for sample_name in samples_with_paralogs_in_greater_than_threshold_genes:
            genes_with_paralogs_handle.write(f'{sample_name}\n')

        genes_with_paralogs_handle.write(f'\nThere are {len(genes_with_paralogs_in_greater_than_threshold_samples)} '
                                         f'genes with paralogs in >= {paralogs_list_threshold_percentage}% of '
                                         f'samples. The gene names are:\n')
        for gene_name in genes_with_paralogs_in_greater_than_threshold_samples:
            genes_with_paralogs_handle.write(f'{gene_name}\n')


def create_paralog_heatmap(paralog_report_filename,
                           figure_length,
                           figure_height,
                           sample_text_size,
                           gene_text_size,
                           heatmap_filename,
                           heatmap_filetype,
                           heatmap_dpi):

    """

    :param paralog_report_filename: filename of the paralog report file *.tsv table
    :param int figure_length: length in inches for heatmap figure
    :param int figure_height: height in inches for heatmap figure
    :param int sample_text_size: size in points for sample name text in heatmap figure
    :param int gene_text_size: size in points for gene name text in heatmap figure
    :param str heatmap_filename: name for the heatmap figure image file
    :param str heatmap_filetype: type of figure to save ('png', 'pdf', 'eps', 'tiff', 'svg')
    :param in heatmap_dpi: Dots per inch (DPI) for the output heatmap image
    :return:
    """

    # Read in the sequence length file:
    df = pd.read_csv(paralog_report_filename, delimiter='\t', )

    # Melt the dataframe to make it suitable for the df.pivot method:
    df = df.melt(id_vars=['Species'], var_name='gene_id', value_name='num_paralogs')

    # Pivot the dataframe for input into the seaborn heatmap function:
    df = df.pivot(index='Species', columns='gene_id', values='num_paralogs')

    # Get figure dimension and label text size based on number of samples and genes:
    fig_length, figure_height, sample_text_size, gene_id_text_size = get_figure_dimensions(df,
                                                                                           figure_length,
                                                                                           figure_height,
                                                                                           sample_text_size,
                                                                                           gene_text_size)

    # Create heatmap:
    sns.set(rc={'figure.figsize': (fig_length, figure_height)})
    sns.set_style('ticks')  # options are: white, dark, whitegrid, darkgrid, ticks
    cmap = 'bone_r'  # sets colour scheme
    heatmap = sns.heatmap(df, vmin=0, cmap=cmap)
    heatmap.tick_params(axis='x', labelsize=gene_id_text_size)
    heatmap.tick_params(axis='y', labelsize=sample_text_size)
    heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0)
    heatmap.set_xlabel("Gene name", fontsize=14, fontweight='bold', labelpad=20)
    heatmap.set_ylabel("Sample name", fontsize=14, fontweight='bold', labelpad=20)
    plt.title("Number of paralog sequences for each gene, for each sample", fontsize=14,
              fontweight='bold', y=1.05)
    # plt.tight_layout()

    # Save heatmap as png file:
    logger.info(f'Saving heatmap as file "{heatmap_filename}.{heatmap_filetype}" at {heatmap_dpi} DPI')
    plt.savefig(f'{heatmap_filename}.{heatmap_filetype}', dpi=heatmap_dpi, bbox_inches='tight')


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
    parser.add_argument('--paralogs_above_threshold_report_filename',
                        help='Specify the filename for the *.txt list of genes with paralogs in '
                             '<paralogs_list_threshold_percentage> number of samples',
                        default='paralogs_above_threshold_report')
    parser.add_argument('--paralogs_list_threshold_percentage',
                        help='Percent of total number of samples and genes that must have paralog warnings to be '
                             'reported in the <genes_with_paralogs.txt> report file. The default of 0.0 means that '
                             'all genes and samples with at least one paralog warning will be reported',
                        type=float,
                        default=0.0)
    parser.add_argument('--heatmap_filename',
                        help='Filename for the output heatmap, saved by default as a *.png file. Defaults to '
                             '"paralog_heatmap"',
                        default='paralog_heatmap')
    parser.add_argument('--figure_length', type=int,
                        help='Length dimension (in inches) for the output heatmap file. Default is '
                             'automatically calculated based on the number of genes',
                        default=None)
    parser.add_argument('--figure_height', type=int,
                        help='Height dimension (in inches) for the output heatmap file. Default is '
                             'automatically calculated based on the number of samples',
                        default=None)
    parser.add_argument('--sample_text_size', type=int,
                        help='Size (in points) for the sample text labels in the output heatmap file. Default is '
                             'automatically calculated based on the number of samples',
                        default=None)
    parser.add_argument('--gene_text_size', type=int,
                        help='Size (in points) for the gene text labels in the output heatmap file. Default is '
                             'automatically calculated based on the number of genes',
                        default=None)
    parser.add_argument('--heatmap_filetype',
                        choices=['png', 'pdf', 'eps', 'tiff', 'svg'],
                        help='File type to save the output heatmap image as. Default is png',
                        default='png')
    parser.add_argument('--heatmap_dpi', type=int,
                        help='Dots per inch (DPI) for the output heatmap image. Default is 300',
                        default='300')

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

    stats_for_tsv_report_all_samples = []
    genes_with_paralogs_sample_lists = []
    for name in namelist:  # iterate over samples
        sample_base_directory_path, sample_directory_name = os.path.split(name)
        if not sample_directory_name:
            sample_base_directory_path, sample_folder_name = os.path.split(sample_base_directory_path)

        stats_for_tsv_report, genes_with_paralogs = retrieve_seqs(sample_base_directory_path,
                                                                  sample_directory_name,
                                                                  target_genes,
                                                                  args.fasta_dir_all,
                                                                  args.fasta_dir_no_chimeras)

        stats_for_tsv_report_all_samples.append(stats_for_tsv_report)
        genes_with_paralogs_sample_lists.append(genes_with_paralogs)

    # Write text statistics report:
    write_paralogs_above_threshold_report(genes_with_paralogs_sample_lists,
                                          args.paralogs_list_threshold_percentage,
                                          namelist,
                                          target_genes,
                                          args.paralogs_above_threshold_report_filename)

    # Write a *.tsv report, i.e. a matrix of sample names vs gene name with sequence counts:
    with open(f'{args.paralog_report_filename}.tsv', 'w') as paralog_report_handle:
        genes = '\t'.join(target_genes)
        paralog_report_handle.write(f'Species\t{genes}\n')
        for stat_list in stats_for_tsv_report_all_samples:
            stats = '\t'.join([str(stat) for stat in stat_list])
            paralog_report_handle.write(f'{stats}\n')

    # Create a heatmap from the *.tsv file:
    create_paralog_heatmap(f'{args.paralog_report_filename}.tsv',
                           args.figure_length,
                           args.figure_height,
                           args.sample_text_size,
                           args.gene_text_size,
                           args.heatmap_filename,
                           args.heatmap_filetype,
                           args.heatmap_dpi)


if __name__ == "__main__":
    standalone()
