#!/usr/bin/env python

"""
This script will retrieve paralog nucleotide (CDS) sequences for a specified gene in all samples located in
namelist.txt. It writes the unaligned sequences to folders, which can be user specified; defaults are 'paralogs_all'
and 'paralogs_no_chimeras'.

If a sample does not have paralogs detected for that gene, the sequence in the FNA directory is retrieved instead.
"""

import os
import sys
import argparse
from Bio import SeqIO
import logging
from collections import defaultdict
import progressbar
import textwrap
import re

from hybpiper.gene_recovery_heatmap import get_figure_dimensions
from hybpiper.retrieve_sequences import get_chimeric_genes_for_sample
from hybpiper import utils
from hybpiper.version import __version__

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
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)

# Setup logger:
logger = logging.getLogger(f'hybpiper.{__name__}')

# Add handlers to the logger
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)  # Default level is 'WARNING'


# Set widget format for progressbar:
widgets = [' ' * 11,
           progressbar.Timer(),
           progressbar.Bar(),
           progressbar.ETA()]


def retrieve_gene_paralogs_from_sample(sampledir_parent,
                                       sample_name,
                                       compressed_sample_dict,
                                       compressed_sample_bool,
                                       gene,
                                       chimera_check_performed):
    """
    Takes a gene name for a given sample, and for each gene produces two *.fasta files:  1) all
    paralog sequences (or *.FNA if no paralogs) from all samples; 2)  all paralog sequences (or non-chimeric *.FNA if
    no paralogs) from all samples; the latter is only written if a chimera check was performed for the sample during
    'hybpiper assemble'.

    :param str sampledir_parent:
    :param str sample_name: name (and expected path relative to cwd) of the sample
    :param dict compressed_sample_dict:
    :param bool compressed_sample_bool:
    :param str gene: name of the gene to recover sequences for
    :param bool chimera_check_performed: True is a chimera check was performed for this sample
    :return list stats_for_report, genes_with_paralogs: lists of sequence counts and gene names for writing reports
    """

    chimeric_genes_to_skip = []
    if chimera_check_performed:
        chimeric_genes_to_skip = get_chimeric_genes_for_sample(sampledir_parent,
                                                               sample_name,
                                                               compressed_sample_dict,
                                                               compressed_sample_bool)

    # Normal recovery of all sequences:
    seqs_to_write_all = []
    seqs_to_write_no_chimeras = []
    num_seqs = 0  # default if no paralogs or *.FNA
    has_paralogs = False

    ####################################################################################################################
    # Get paths to paralogs file (might not be present if no paralogs) and *.FNA file; sample dir as root:
    ####################################################################################################################
    paralog_fasta_file_path = os.path.join(sample_name,
                                           gene,
                                           sample_name,
                                           'paralogs',
                                           f'{gene}_paralogs.fasta')

    fna_fasta_file_path = os.path.join(sample_name,
                                       gene,
                                       sample_name,
                                       'sequences',
                                       'FNA',
                                       f'{gene}.FNA')

    # Full paths:
    full_paralog_fasta_file_path = f'{sampledir_parent}/{paralog_fasta_file_path}'
    full_fna_fasta_file_path = f'{sampledir_parent}/{fna_fasta_file_path}'


    ####################################################################################################################
    # Recover paralog sequences if present:
    ####################################################################################################################
    if compressed_sample_bool:
        if paralog_fasta_file_path in compressed_sample_dict[sample_name]:
            paralog_fasta_file_path_size = compressed_sample_dict[sample_name][paralog_fasta_file_path]

            if paralog_fasta_file_path_size  == 0:
                logger.warning(f'{"[WARNING]:":10} File {paralog_fasta_file_path} exists, but is empty!')
            else:
                seqs_to_write = utils.get_compressed_seqrecords(sample_name,
                                                                sampledir_parent,
                                                                paralog_fasta_file_path)
                seqs_to_write_all.extend(seqs_to_write)
                num_seqs = len(seqs_to_write)
                has_paralogs = True
    else:
        if os.path.isfile(full_paralog_fasta_file_path):
            if os.path.getsize(full_paralog_fasta_file_path) == 0:
                logger.warning(f'{"[WARNING]:":10} File {paralog_fasta_file_path} exists, but is empty!\n')
            else:
                seqs_to_write = [x for x in SeqIO.parse(full_paralog_fasta_file_path, 'fasta')]
                seqs_to_write_all.extend(seqs_to_write)
                num_seqs = len(seqs_to_write)
                has_paralogs = True

    ####################################################################################################################
    # If there are no paralogs, recover the nucleotide *.FNA sequence if present:
    ####################################################################################################################
    if not has_paralogs:
        if compressed_sample_bool:
            if fna_fasta_file_path in compressed_sample_dict[sample_name]:
                fna_fasta_file_path_size = compressed_sample_dict[sample_name][fna_fasta_file_path]

                if fna_fasta_file_path_size == 0:
                    logger.warning(f'{"[WARNING]:":10} File {fna_fasta_file_path} exists, but is empty!')

                else:
                    seqs_to_write = utils.get_compressed_seqrecords(sample_name,
                                                                    sampledir_parent,
                                                                    fna_fasta_file_path)
                    seqs_to_write_all.extend(seqs_to_write)
                    num_seqs = 1
        else:
            if os.path.isfile(full_fna_fasta_file_path):
                if os.path.getsize(full_fna_fasta_file_path) == 0:
                    logger.warning(f'{"[WARNING]:":10} File {fna_fasta_file_path} exists, but is empty!\n')
                else:
                    seqs_to_write = SeqIO.read(full_fna_fasta_file_path, 'fasta')
                    seqs_to_write_all.append(seqs_to_write)
                    num_seqs = 1

    # Now skip any putative chimeric stitched contig sequences (if chimera check was performed for this sample). Note
    # that if paralogs are present they are recovered regardless of whether the gene is in chimeric_genes_to_skip,
    # as paralog sequences are derived from single Exonerate hits only:
    if chimera_check_performed and gene in chimeric_genes_to_skip:
        if compressed_sample_bool:
            if paralog_fasta_file_path in compressed_sample_dict[sample_name]:
                paralog_fasta_file_path_size = compressed_sample_dict[sample_name][paralog_fasta_file_path]

                if paralog_fasta_file_path_size == 0:
                    logger.warning(f'{"[WARNING]:":10} File {paralog_fasta_file_path} exists, but is empty!')
                else:
                    seqs_to_write = utils.get_compressed_seqrecords(sample_name,
                                                                    sampledir_parent,
                                                                    paralog_fasta_file_path)
                    seqs_to_write_no_chimeras.extend(seqs_to_write)
            else:
                pass  # don't get the *.FNA sequence for this sample/gene combination

        else:
            if os.path.isfile(full_paralog_fasta_file_path):
                if os.path.getsize(full_paralog_fasta_file_path) == 0:
                    logger.warning(f'{"[WARNING]:":10} File {paralog_fasta_file_path} exists, but is empty!\n')
                else:
                    seqs_to_write_no_chimeras.extend([x for x in SeqIO.parse(full_paralog_fasta_file_path,
                                                                             'fasta')])
            elif os.path.isfile(full_fna_fasta_file_path):
                if os.path.getsize(full_fna_fasta_file_path) == 0:
                    logger.warning(f'{"[WARNING]:":10} File {fna_fasta_file_path} exists, but is empty!\n')
                else:
                    seqs_to_write_no_chimeras.append(SeqIO.read(full_fna_fasta_file_path, 'fasta'))

    return (num_seqs,
            has_paralogs,
            seqs_to_write_all,
            seqs_to_write_no_chimeras)


def write_paralogs_above_threshold_report(gene_with_paralogs_to_sample_list_dict, paralogs_list_threshold_percentage,
                                          namelist, target_genes, paralogs_above_threshold_report_filename):
    """
    Writes a *.txt report list gene names and sample names for genes and samples that have paralogs in >=

    :param dict gene_with_paralogs_to_sample_list_dict: dict of lists e.g. {'gene006': ['EG98'], 'gene026': ['EG30',
    'EG98'], etc. }
    :param float paralogs_list_threshold_percentage: percentage threshold used for report calculations
    :param list namelist: list of sample names
    :param list target_genes: list of target gene names
    :param str paralogs_above_threshold_report_filename: file name for the report file
    :return:
    """

    # Calculate number of genes with warnings in >= paralogs_list_threshold_percentage samples:
    genes_with_paralogs_in_greater_than_threshold_samples = []

    for gene, sample_name_list in gene_with_paralogs_to_sample_list_dict.items():
        if len(sample_name_list) / len(namelist) >= paralogs_list_threshold_percentage:  # Check percentage of samples
            genes_with_paralogs_in_greater_than_threshold_samples.append(gene)

    logger.info(f'{"[INFO]:":10} There are {len(genes_with_paralogs_in_greater_than_threshold_samples)} genes with '
                f'paralogs in >= {paralogs_list_threshold_percentage}% of samples. The gene names have been written '
                f'to the file "{paralogs_above_threshold_report_filename}.txt"')

    # Calculate number of samples with warnings in >= paralogs_list_threshold_percentage genes:
    samples_with_paralogs_in_greater_than_threshold_genes = []
    sample_to_paralogous_genes_count_dict = defaultdict(int)

    for gene, sample_name_list in gene_with_paralogs_to_sample_list_dict.items():
        for sample_name in sample_name_list:
            sample_to_paralogous_genes_count_dict[sample_name] += 1

    for sample, paralogous_genes_count in sample_to_paralogous_genes_count_dict.items():
        if paralogous_genes_count / len(target_genes) >= paralogs_list_threshold_percentage:
            samples_with_paralogs_in_greater_than_threshold_genes.append(sample)

    logger.info(f'{"[INFO]:":10} There are {len(samples_with_paralogs_in_greater_than_threshold_genes)} samples with '
                f'paralogs in >= {paralogs_list_threshold_percentage}% of genes. The sample names have been written '
                f'to the file "{paralogs_above_threshold_report_filename}.txt"')

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
                           heatmap_dpi,
                           no_xlabels,
                           no_ylabels):

    """
    Creates a heatmap showing the number of sequences recovered for each gene, for each sample.

    :param paralog_report_filename: filename of the paralog report file *.tsv table
    :param int figure_length: length in inches for heatmap figure
    :param int figure_height: height in inches for heatmap figure
    :param int sample_text_size: size in points for sample name text in heatmap figure
    :param int gene_text_size: size in points for gene name text in heatmap figure
    :param str heatmap_filename: name for the heatmap figure image file
    :param str heatmap_filetype: type of figure to save ('png', 'pdf', 'eps', 'tiff', 'svg')
    :param int heatmap_dpi: Dots per inch (DPI) for the output heatmap image
    :param bool no_xlabels: if True, don't draw x-axis labels on saved figure
    :param bool no_ylabels: if True, don't draw y-axis labels on saved figure
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

    # Check that figure won't be greater than the maximum pixels allowed (65536) in either dimension, and resize to
    # 400 inches / 100 DPI if it is. Note that even if the dimensions are less than 65536, a large dataset can still
    #  on render partially (https://stackoverflow.com/questions/64393779/how-to-render-a-heatmap-for-a-large-array):

    figure_length_pixels = fig_length * heatmap_dpi
    figure_height_pixels = figure_height * heatmap_dpi

    if figure_length_pixels >= 65536:

        fig_length = 400
        heatmap_dpi = 100

        fill = textwrap.fill(
            f'{"[INFO]:":10} The large number of loci in this analysis means that the auto-calculated figure length '
            f'({fig_length:.2f} inches / {figure_length_pixels} pixels) is larger than the maximum allowed size '
            f'(65536 pixels). Figure length has been set to 400 inches, and DPI has been set to 100 '
            f'({400 * heatmap_dpi:.2f} pixels). If you find that the heatmap is only partially rendered in the '
            f'saved file, try reducing the DPI further via the --heatmap_dpi parameter, and/or the figure length via '
            f'the --figure_length parameter.',
            width=90, subsequent_indent=' ' * 11)

        logger.info(fill)

    if figure_height_pixels >= 65536:

        figure_height = 400
        heatmap_dpi = 100

        fill = textwrap.fill(
            f'{"[INFO]:":10} The large number of samples in this analysis means that the auto-calculated figure height '
            f'({figure_height:.2f} inches / {figure_height_pixels} pixels) is larger than the maximum allowed size '
            f'(65536 pixels). Figure height has been set to 400 inches, and DPI has been set to 100 '
            f'({400 * heatmap_dpi:.2f} pixels). If you find that the heatmap is only partially rendered in the '
            f'saved file, try reducing the DPI further via the --heatmap_dpi parameter, and/or the figure length via '
            f'the --figure_length parameter.',
            width=90, subsequent_indent=' ' * 11)

        logger.info(fill)

    # Create heatmap:
    sns.set(rc={'figure.figsize': (fig_length, figure_height)})
    sns.set_style('ticks')  # options are: white, dark, whitegrid, darkgrid, ticks
    cmap = 'bone_r'  # sets colour scheme
    heatmap = sns.heatmap(df, vmin=0, cmap=cmap, xticklabels=1, yticklabels=1,
                          cbar_kws={"orientation": "vertical", "pad": 0.01})
    heatmap.tick_params(axis='x', labelsize=gene_id_text_size)
    heatmap.tick_params(axis='y', labelsize=sample_text_size)
    heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0)
    heatmap.set_xlabel("Gene name", fontsize=14, fontweight='bold', labelpad=20)
    heatmap.set_ylabel("Sample name", fontsize=14, fontweight='bold', labelpad=20)
    plt.title("Number of paralog sequences for each gene, for each sample", fontsize=14,
              fontweight='bold', y=1.05)

    # Remove x-axis and y-axis labels if flags provided:
    if no_xlabels:
        heatmap.set(xticks=[])

    if no_ylabels:
        heatmap.set(yticks=[])

    # Save heatmap as png file:
    logger.info(f'{"[INFO]:":10} Saving heatmap as file "{heatmap_filename}.{heatmap_filetype}" at {heatmap_dpi} DPI')
    plt.savefig(f'{heatmap_filename}.{heatmap_filetype}', dpi=heatmap_dpi, bbox_inches='tight')


def standalone():
    """
    Used when this module is run as a stand-alone script. Parses command line arguments and runs function main().:
    """

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('namelist',
                        help='Text file containing list of HybPiper output directories, one per line.')

    group_1 = parser.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna', dest='targetfile_dna',
                         help='FASTA file containing DNA target sequences for each gene. Used to extract unique gene '
                              'names for paralog recovery. If there are multiple targets for a gene, the id must be '
                              'of the form: >Taxon-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa', dest='targetfile_aa',
                         help='FASTA file containing amino-acid target sequences for each gene. Used to extract '
                              'unique gene names for paralog recovery. If there are multiple targets for a gene, '
                              'the id must be of the form: >Taxon-geneName')
    parser.add_argument('--hybpiper_dir', help='Specify directory containing HybPiper output',
                        default=None)
    parser.add_argument('--fasta_dir_all',
                        help='Specify directory for output FASTA files (ALL). Default is "paralogs_all".',
                        default='paralogs_all')
    parser.add_argument('--fasta_dir_no_chimeras',
                        help='Specify directory for output FASTA files (no putative chimeric sequences). Default is '
                             '"fasta_dir_no_chimeras".',
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
    parser.add_argument('--no_xlabels',
                        action='store_true',
                        default=False,
                        help='If supplied, do not render labels for x-axis (loci) in the saved heatmap figure')
    parser.add_argument('--no_ylabels',
                        action='store_true',
                        default=False,
                        help='If supplied, do not render labels for y-axis (samples) in the saved heatmap figure')

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

    logger.info(f'{"[INFO]:":10} Recovering paralog sequences from the HybPiper run(s)...')

    # Set target file name:
    targetfile = args.targetfile_dna if args.targetfile_dna else args.targetfile_aa

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

    ####################################################################################################################
    # Parse namelist and check for the presence of the corresponding sample directories or *.tar.gz files:
    ####################################################################################################################
    list_of_sample_names = []

    with open(args.namelist, 'r') as namelist_handle:
        for line in namelist_handle.readlines():
            sample_name = line.rstrip()
            if sample_name:
                list_of_sample_names.append(sample_name)
                if re.search('/', sample_name):
                    sys.exit(f'{"[ERROR]:":10} A sample name must not contain '
                             f'forward slashes. The file {args.namelist} contains: {sample_name}')

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

        fill = utils.fill_forward_slash(f'{"[WARNING]:":10} File {args.namelist} contains samples not found in '
                                        f'directory "{sampledir_parent}". The missing samples are:',
                                        width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                        break_on_forward_slash=True)

        logger.warning(f'{fill}\n')

        for name in samples_missing:
            logger.warning(f'{" " * 10} {name}')
        logger.warning('')

    list_of_sample_names = [x for x in list_of_sample_names if x not in samples_missing]

    ########################################################################################################
    # Determine whether a chimera check was performed for each sample during 'hybpiper assemble':
    ########################################################################################################
    chimera_check_samples_dict = dict()

    for sample_name in list_of_sample_names:

        # Check if the sample directory is a compressed tarball:
        compressed_sample_bool = True if sample_name in compressed_sample_dict else False

        chimera_check_performed_file = f'{sample_name}/{sample_name}_chimera_check_performed.txt'  # sample dir as root

        chimera_file_fill = textwrap.fill(f'{"[ERROR]:":10} No file "{chimera_check_performed_file}" '
                                          f'found. If you are running `hybpiper paralog_retriever` from Hybpiper '
                                          f'version >=2.2.0 on a sample that was assembled using HybPiper version '
                                          f'<2.2.0, please create a text file in the main sample directory for sample '
                                          f'{sample_name} named "{chimera_check_performed_file}", containing '
                                          f'the text "True" (no quotation marks), and run `hybpiper paralog_retriever` '
                                          f'again.',
                                          width=90, subsequent_indent=" " * 11)

        if compressed_sample_bool:
            try:
                compressed_sample_bool_lines = utils.get_compressed_file_lines(sample_name,
                                                                               sampledir_parent,
                                                                               chimera_check_performed_file)
            except KeyError:
                logger.error(chimera_file_fill)
                sys.exit(1)

            chimera_check_bool = compressed_sample_bool_lines[0]

        else:
            try:
                with open(f'{sampledir_parent}/{chimera_check_performed_file}', 'r') as chimera_check_handle:
                    chimera_check_bool = chimera_check_handle.read().rstrip()
            except FileNotFoundError:
                logger.error(chimera_file_fill)
                sys.exit(1)

        if chimera_check_bool == 'True':
            chimera_check_samples_dict[sample_name] = True
        elif chimera_check_bool == 'False':
            chimera_check_samples_dict[sample_name] = False
        else:
            raise ValueError(f'chimera_check_bool is: {chimera_check_bool} for sample {sample_name}')

    # Check is at least one sample had a chimera check performed:
    chimera_check_performed_for_at_least_one_sample = False
    for sample in chimera_check_samples_dict:
        if chimera_check_samples_dict[sample]:
            chimera_check_performed_for_at_least_one_sample = True
            break

    ####################################################################################################################
    # Set output directories:
    ####################################################################################################################
    if chimera_check_performed_for_at_least_one_sample:
        logger.info(f'{"[INFO]:":10} Creating directory: {args.fasta_dir_no_chimeras}')
        if not os.path.isdir(args.fasta_dir_no_chimeras):
            os.mkdir(args.fasta_dir_no_chimeras)

    logger.info(f'{"[INFO]:":10} Creating directory: {args.fasta_dir_all}')
    if not os.path.isdir(args.fasta_dir_all):
        os.mkdir(args.fasta_dir_all)

    ####################################################################################################################
    # Get gene names parsed from a target file.
    ####################################################################################################################
    target_genes = sorted(list(set([x.id.split('-')[-1] for x in SeqIO.parse(targetfile, 'fasta')])))

    logger.info(f'{"[INFO]:":10} Searching for paralogs for {len(list_of_sample_names)} samples, '
                f'{len(target_genes)} genes...')

    ####################################################################################################################
    # Warn user if some samples didn't have a chimera check performed during `hybpiper assemble`:
    ####################################################################################################################
    samples_with_no_chimera_check_performed = [sample for sample in chimera_check_samples_dict if not
                                               chimera_check_samples_dict[sample]]

    if chimera_check_performed_for_at_least_one_sample:  # else don't warn - chimera-filtered output folder not present
        fill = textwrap.fill(f'{"[WARNING]:":10} A chimera check was not performed during "hybpiper assemble" for '
                             f'the following samples:',
                             width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
        logger.warning(f'\n{fill}\n')

        for sample in samples_with_no_chimera_check_performed:
            logger.warning(f'{" " * 10} {sample}')

        logger.warning(f'\n{" " * 10} No putative chimeric sequences will be skipped for these samples!\n')

    ####################################################################################################################
    # Create dictionaries to capture data for reports/heatmap:
    ####################################################################################################################
    sample_to_gene_paralog_count_dict = defaultdict(dict)
    gene_with_paralogs_to_sample_list_dict = defaultdict(list)

    ####################################################################################################################
    # Iterate over each gene and capture data for all samples:
    ####################################################################################################################
    for gene in progressbar.progressbar(target_genes, max_value=len(target_genes), min_poll_interval=30,
                                        widgets=widgets):
        sequences_to_write_all = []
        sequences_to_write_no_chimeras = []

        ################################################################################################################
        # Iterate over each sample:
        ################################################################################################################
        for sample_name in list_of_sample_names:

            # Check if the sample directory is a compressed tarball:
            compressed_sample_bool = False
            if sample_name in compressed_sample_dict:
                compressed_sample_bool = True

            chimera_check_performed = chimera_check_samples_dict[sample_name]

            (num_seqs,
             has_paralogs,
             seqs_to_write_all,
             seqs_to_write_no_chimeras) = retrieve_gene_paralogs_from_sample(sampledir_parent,
                                                                             sample_name,
                                                                             compressed_sample_dict,
                                                                             compressed_sample_bool,
                                                                             gene,
                                                                             chimera_check_performed)

            # Add the number of sequences for this gene to the sample in the sample_to_gene_paralog_count_dict:
            sample_to_gene_paralog_count_dict[sample_name][gene] = num_seqs

            # If the gene has paralogs for this sample, add the sample name to the list for the gene in
            # gene_with_paralogs_to_sample_list_dict:
            if has_paralogs:
                gene_with_paralogs_to_sample_list_dict[gene].append(sample_name)

            # Add any paralog or *.FNA sequences to the write-to-file lists for this gene:
            sequences_to_write_all.extend(seqs_to_write_all)
            sequences_to_write_no_chimeras.extend(seqs_to_write_no_chimeras)

        # Write the 'all' and 'no chimeras' fasta files for this gene:
        with open(f'{args.fasta_dir_all}/{gene}_paralogs_all.fasta', 'w') as all_seqs_handle:
            SeqIO.write(sequences_to_write_all, all_seqs_handle, 'fasta')

        if chimera_check_performed_for_at_least_one_sample:  # chimera-filtered output folder will be present
            with open(f'{args.fasta_dir_no_chimeras}/{gene}_paralogs_no_chimeras.fasta', 'w') as no_chimera_seqs_handle:
                SeqIO.write(sequences_to_write_no_chimeras, no_chimera_seqs_handle, 'fasta')

    # Write a *.tsv report, i.e. a matrix of sample names vs gene name with sequence counts:
    with open(f'{args.paralog_report_filename}.tsv', 'w') as paralog_report_handle:
        genes = '\t'.join(target_genes)
        paralog_report_handle.write(f'Species\t{genes}\n')
        for sample, gene_paralog_count_dict in sample_to_gene_paralog_count_dict.items():
            gene_counts = '\t'.join([str(count) for count in gene_paralog_count_dict.values()])
            paralog_report_handle.write(f'{sample}\t{gene_counts}\n')

    # Create a heatmap from the *.tsv file:
    if args.no_heatmap:
        logger.info(f'{"[INFO]:":10} Option "--no_heatmap" provided. No paralog heatmap will be created!')
    else:
        logger.info(f'{"[INFO]:":10} Creating paralog heatmap...')
        create_paralog_heatmap(f'{args.paralog_report_filename}.tsv',
                               args.figure_length,
                               args.figure_height,
                               args.sample_text_size,
                               args.gene_text_size,
                               args.heatmap_filename,
                               args.heatmap_filetype,
                               args.heatmap_dpi,
                               args.no_xlabels,
                               args.no_ylabels)

    # Write text statistics report:
    write_paralogs_above_threshold_report(gene_with_paralogs_to_sample_list_dict,
                                          args.paralogs_list_threshold_percentage,
                                          list_of_sample_names,
                                          target_genes,
                                          args.paralogs_above_threshold_report_filename)


if __name__ == "__main__":
    standalone()
