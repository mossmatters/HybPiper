#!/usr/bin/env python

"""
Takes the seq_lengths.txt file (output from running 'hybpiper get_seq_lengths') as input.

For each sample and for each gene, calculates the percentage length recovered. This percentage is calculated as a
fraction of the mean length for representative gene sequences in the target file provided.

Generates a heatmap of percentage length recovery for each sample and each gene.

"""

import sys
import argparse
import os
import logging
import textwrap

from hybpiper import utils
from hybpiper.version import __version__

# Import non-standard-library modules:

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

########################################################################################################################
########################################################################################################################
# Get current working directory and host name:

cwd = os.getcwd()


########################################################################################################################
########################################################################################################################
# Define functions:

def get_figure_dimensions(df, figure_length, figure_height, sample_text_size, gene_text_size):
    """
    Takes a dataframe and returns figure length (inches), figure height (inches), sample_text_size and gene_id_text_size
    values based on the number of samples and genes in the seq_lengths.txt input file provided.

    :param pandas.core.frame.DataFrame df: pandas dataframe of seq_lengths.txt after filtering and pivot
    :param NoneType or int figure_length: if provided, dimension (in inches) for the figure length
    :param NoneType or int figure_height: if provided, dimension (in inches) for the figure height
    :param NoneType or int sample_text_size: if provided, dimension (in inches) for the figure sample text size
    :param NoneType or int gene_text_size: if provided, dimension (in inches) for the figure gene_id text size
    :return float fig_length, figure_height, sample_text_size, gene_id_text_size:
    """

    num_samples = len(df.index)
    num_genes = len(df.columns)

    logger.info(f'{"[INFO]:":10} Number of samples in input lengths file is: {num_samples}')
    logger.info(f'{"[INFO]:":10} Number of genes in input lengths file is: {num_genes}')

    # Set default text label size (in points) unless specified at the command line:
    sample_text_size = sample_text_size if sample_text_size else 10
    gene_id_text_size = gene_text_size if gene_text_size else 10

    # Set figure height dimensions for a given number of samples:
    figure_height = figure_height if figure_height else num_samples / 3

    # Set figure length dimensions for a given number of genes (i.e. number of unique genes in target file):
    fig_length = figure_length if figure_length else num_genes / 3

    logger.info(f'{"[INFO]:":10} figure_length: {fig_length:.2f} inches, figure_height: {figure_height:.2f} inches, '
                f'sample_text_size: {sample_text_size} points, gene_id_text_size: {gene_id_text_size} points')

    return fig_length, figure_height, sample_text_size, gene_id_text_size


########################################################################################################################
########################################################################################################################
# Run script:

def standalone():
    """
    Used when this module is run as a stand-alone script. Parses command line arguments and runs function main().:
    """

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('seq_lengths_file',
                        help="Filename for the seq_lengths file (output by the 'hybpiper stats' command)")
    parser.add_argument('--heatmap_filename',
                        default='recovery_heatmap',
                        help='Filename for the output heatmap, saved by default as a *.png file. Defaults to '
                             '"recovery_heatmap"')
    parser.add_argument('--figure_length',
                        type=int,
                        default=None,
                        help='Length dimension (in inches) for the output heatmap file. Default is '
                             'automatically calculated based on the number of genes')
    parser.add_argument('--figure_height',
                        type=int,
                        default=None,
                        help='Height dimension (in inches) for the output heatmap file. Default is '
                             'automatically calculated based on the number of samples')
    parser.add_argument('--sample_text_size',
                        type=int,
                        default=None,
                        help='Size (in points) for the sample text labels in the output heatmap file. Default is '
                             'automatically calculated based on the number of samples')
    parser.add_argument('--gene_text_size',
                        type=int,
                        default=None,
                        help='Size (in points) for the gene text labels in the output heatmap file. Default is '
                             'automatically calculated based on the number of genes')
    parser.add_argument('--heatmap_filetype', choices=['png', 'pdf', 'eps', 'tiff', 'svg'],
                        help='File type to save the output heatmap image as. Default is png',
                        default='png')
    parser.add_argument('--heatmap_dpi',
                        type=int,
                        default=100,
                        help='Dots per inch (DPI) for the output heatmap image. Default is %(default)d')
    parser.add_argument('--no_xlabels',
                        action='store_true',
                        default=False,
                        help='If supplied, do not render labels for x-axis (loci) in the saved heatmap figure')
    parser.add_argument('--no_ylabels',
                        action='store_true',
                        default=False,
                        help='If supplied, do not render labels for y-axis (samples) in the saved heatmap figure')

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

    logger.info(f'{"[INFO]:":10} Creating a gene-recovery heatmap for the HybPiper run(s)...')

    ####################################################################################################################
    # Check for presence of required input files:
    ####################################################################################################################
    logger.info(f'{"[INFO]:":10} The following file of sequence lengths was provided: "{args.seq_lengths_file}".')
    if not utils.file_exists_and_not_empty(args.seq_lengths_file):
        sys.exit(f'{"[ERROR]:":10} File {args.seq_lengths_file} is missing or empty, exiting!')

    ####################################################################################################################
    # Read in the sequence length file and process:
    ####################################################################################################################
    df = pd.read_csv(args.seq_lengths_file, delimiter='\t', )

    df = df.astype('object')  # https://pandas.pydata.org/docs/whatsnew/v2.1.0.html#deprecated-silent-upcasting-in-setitem-like-series-operations

    # For each sample, divide each gene length by the MeanLength value for that gene:
    df.loc[:, df.columns[1]:] = df.loc[:, df.columns[1]:].div(df.iloc[0][df.columns[1]:])

    # For each length ratio, if the value is greater than 1, assign it to 1:
    df.where(df.loc[:, df.columns[1]:] < 1, 1, inplace=True)

    # Drop the gene MeanLengths row:
    df.drop(labels=0, axis=0, inplace=True)

    # Melt the dataframe to make it suitable for the df.pivot method:
    df = df.melt(id_vars=['Species'], var_name='gene_id', value_name='percentage_recovery')

    # Change percentage values to numeric:
    df["percentage_recovery"] = df["percentage_recovery"].apply(pd.to_numeric)

    # Pivot the dataframe for input into the seaborn heatmap function:
    df = df.pivot(index='Species', columns='gene_id', values='percentage_recovery')

    ####################################################################################################################
    # Get figure dimension and label text size based on number of samples and genes:
    ####################################################################################################################
    fig_length, figure_height, sample_text_size, gene_id_text_size = get_figure_dimensions(df,
                                                                                           args.figure_length,
                                                                                           args.figure_height,
                                                                                           args.sample_text_size,
                                                                                           args.gene_text_size)

    # Check that figure won't be greater than the maximum pixels allowed (65536) in either dimension, and resize to
    # 400 inches / 100 DPI if it is. Note that even if the dimensions are less than 65536, a large dataset can still
    #  on render partially (https://stackoverflow.com/questions/64393779/how-to-render-a-heatmap-for-a-large-array):

    figure_length_pixels = fig_length * args.heatmap_dpi
    figure_height_pixels = figure_height * args.heatmap_dpi

    if figure_length_pixels >= 65536:

        fig_length = 400
        args.heatmap_dpi = 100

        fill = textwrap.fill(
            f'{"[INFO]:":10} The large number of loci in this analysis means that the auto-calculated figure length '
            f'({fig_length:.2f} inches / {figure_length_pixels} pixels) is larger than the maximum allowed size '
            f'(65536 pixels). Figure length has been set to 400 inches, and DPI has been set to 100 '
            f'({400 * args.heatmap_dpi:.2f} pixels). If you find that the heatmap is only partially rendered in the '
            f'saved file, try reducing the DPI further via the --heatmap_dpi parameter, and/or the figure length via '
            f'the --figure_length parameter.',
            width=90, subsequent_indent=' ' * 11)

        logger.info(fill)

    if figure_height_pixels >= 65536:

        figure_height = 400
        args.heatmap_dpi = 100

        fill = textwrap.fill(
            f'{"[INFO]:":10} The large number of samples in this analysis means that the auto-calculated figure height '
            f'({figure_height:.2f} inches / {figure_height_pixels} pixels) is larger than the maximum allowed size '
            f'(65536 pixels). Figure height has been set to 400 inches, and DPI has been set to 100 '
            f'({400 * args.heatmap_dpi:.2f} pixels). If you find that the heatmap is only partially rendered in the '
            f'saved file, try reducing the DPI further via the --heatmap_dpi parameter, and/or the figure length via '
            f'the --figure_length parameter.',
            width=90, subsequent_indent=' ' * 11)

        logger.info(fill)

    ####################################################################################################################
    # Create heatmap:
    ####################################################################################################################
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
    plt.title("Percentage length recovery for each gene, relative to mean of targetfile references", fontsize=14,
              fontweight='bold', y=1.05)

    # Remove x-axis and y-axis labels if flags provided:
    if args.no_xlabels:
        heatmap.set(xticks=[])

    if args.no_ylabels:
        heatmap.set(yticks=[])

    ####################################################################################################################
    # Save heatmap to file:
    ####################################################################################################################
    logger.info(f'{"[INFO]:":10} Saving heatmap as file "{args.heatmap_filename}.{args.heatmap_filetype}" at'
                f' {args.heatmap_dpi} DPI')
    plt.savefig(f'{args.heatmap_filename}.{args.heatmap_filetype}', dpi=args.heatmap_dpi, bbox_inches='tight')


########################################################################################################################
########################################################################################################################

if __name__ == '__main__':
    standalone()

########################################################################################################################
########################################################################################################################
