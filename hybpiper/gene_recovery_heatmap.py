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
console_handler = logging.StreamHandler(sys.stderr)
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
    figure_height = figure_height if figure_height else num_samples/4

    # Set figure length dimensions for a given number of genes (i.e. number of unique genes in target file):
    fig_length = figure_length if figure_length else num_genes/4

    logger.info(f'{"[INFO]:":10} figure_length: {fig_length} inches, figure_height: {figure_height} inches, '
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
                        help='Filename for the output heatmap, saved by default as a *.png file. Defaults to '
                             '"recovery_heatmap"',
                        default='recovery_heatmap')
    parser.add_argument('--figure_length', type=int,
                        help='Length dimension (in inches) for the output heatmap file. Default is '
                             'automatically calculated based on the number of genes', default=None)
    parser.add_argument('--figure_height', type=int,
                        help='Height dimension (in inches) for the output heatmap file. Default is '
                             'automatically calculated based on the number of samples', default=None)
    parser.add_argument('--sample_text_size', type=int,
                        help='Size (in points) for the sample text labels in the output heatmap file. Default is '
                             'automatically calculated based on the number of samples', default=None)
    parser.add_argument('--gene_text_size', type=int,
                        help='Size (in points) for the gene text labels in the output heatmap file. Default is '
                             'automatically calculated based on the number of genes', default=None)
    parser.add_argument('--heatmap_filetype', choices=['png', 'pdf', 'eps', 'tiff', 'svg'],
                        help='File type to save the output heatmap image as. Default is png', default='png')
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

    # logger.info(f'Running {__name__} with: {args}')

    if args.seq_lengths_file and not os.path.exists(args.seq_lengths_file):
        logger.info(f'Can not find file "{args.seq_lengths_file}". Is it in the current working directory?')

    # Read in the sequence length file:
    df = pd.read_csv(args.seq_lengths_file, delimiter='\t', )

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

    # Get figure dimension and label text size based on number of samples and genes:
    fig_length, figure_height, sample_text_size, gene_id_text_size = get_figure_dimensions(df,
                                                                                           args.figure_length,
                                                                                           args.figure_height,
                                                                                           args.sample_text_size,
                                                                                           args.gene_text_size)

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
    plt.title("Percentage length recovery for each gene, relative to mean of targetfile references", fontsize=14,
              fontweight='bold', y=1.05)
    # plt.tight_layout()

    # Save heatmap as png file:
    logger.info(f'{"[INFO]:":10} Saving heatmap as file "{args.heatmap_filename}.{args.heatmap_filetype}" at'
                f' {args.heatmap_dpi} DPI')
    plt.savefig(f'{args.heatmap_filename}.{args.heatmap_filetype}', dpi=args.heatmap_dpi, bbox_inches='tight')


########################################################################################################################
########################################################################################################################

if __name__ == '__main__':
    standalone()

########################################################################################################################
########################################################################################################################
