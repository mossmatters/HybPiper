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
import multiprocessing
from multiprocessing import Manager
from concurrent.futures.process import ProcessPoolExecutor
from concurrent.futures import wait, as_completed
import traceback
import tarfile

from hybpiper.gene_recovery_heatmap import get_figure_dimensions
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


def write_paralogs_above_threshold_report(gene_with_paralogs_to_sample_list_dict,
                                          paralogs_list_threshold_percentage,
                                          namelist,
                                          target_genes,
                                          paralogs_above_threshold_report_filename):
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


def parse_sample(sample_name,
                 sampledir_parent,
                 compressed_samples_set,
                 target_genes,
                 lock,
                 counter):
    """

    :param str sample_name: name of the sample (no ".tar.gz")
    :param str sampledir_parent: path of the parent directory containing HybPiper output
    :param set compressed_samples_set: set of sample names that are compressed (i.e. *.tar.gz)
    :param list target_genes: list of unique gene names in the target file
    :param multiprocessing.managers.AcquirerProxy lock:
    :param multiprocessing.managers.ValueProxy counter:
    :return:
    """

    compressed_bool = True if sample_name in compressed_samples_set else False
    warning_messages_list = []
    sample_dict = dict()
    sample_dict['paralog_seqrecords'] = dict()
    sample_dict['fna_seqrecords'] = dict()

    ####################################################################################################################
    # Set paths of files to get/process for the sample:
    ####################################################################################################################
    locus_paralog_paths = []
    locus_fna_paths = []

    for gene_name in target_genes:
        paralog_fasta_file_path = os.path.join(sample_name,
                                               gene_name,
                                               sample_name,
                                               'paralogs',
                                               f'{gene_name}_paralogs.fasta')

        locus_paralog_paths.append(paralog_fasta_file_path)

        fna_fasta_file_path = os.path.join(sample_name,
                                           gene_name,
                                           sample_name,
                                           'sequences',
                                           'FNA',
                                           f'{gene_name}.FNA')

        locus_fna_paths.append(fna_fasta_file_path)

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
                    if tarinfo.name in locus_paralog_paths:
                        gene_name = os.path.basename(tarinfo.name).split('_paralogs')[0]
                        sample_dict['paralog_seqrecords'][gene_name] = []

                        seqrecords = utils.get_compressed_seqrecords(tarfile_handle,
                                                                     tarinfo)
                        for seqrecord in seqrecords:
                            sample_dict['paralog_seqrecords'][gene_name].append(seqrecord)

                    if tarinfo.name in locus_fna_paths:
                        gene_name = os.path.basename(tarinfo.name).split('.')[0]
                        sample_dict['fna_seqrecords'][gene_name] = []

                        seqrecords = utils.get_compressed_seqrecords(tarfile_handle,
                                                                     tarinfo)
                        for seqrecord in seqrecords:
                            sample_dict['fna_seqrecords'][gene_name].append(seqrecord)

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
        for locus_paralog_path in locus_paralog_paths:
            locus_path_full = f'{sampledir_parent}/{locus_paralog_path}'
            gene_name = os.path.split(locus_paralog_path)[1].split('_paralog')[0]

            if utils.file_exists_and_not_empty(locus_path_full):
                sample_dict['paralog_seqrecords'][gene_name] = []
                seqrecords = SeqIO.parse(locus_path_full, 'fasta')

                for seqrecord in seqrecords:
                    sample_dict['paralog_seqrecords'][gene_name].append(seqrecord)

        for locus_fna_path in locus_fna_paths:
            locus_path_full = f'{sampledir_parent}/{locus_fna_path}'
            gene_name = os.path.split(locus_fna_path)[1].split('.')[0]

            if utils.file_exists_and_not_empty(locus_path_full):
                sample_dict['fna_seqrecords'][gene_name] = []
                seqrecords = SeqIO.parse(locus_path_full, 'fasta')

                for seqrecord in seqrecords:
                    sample_dict['fna_seqrecords'][gene_name].append(seqrecord)

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

        if utils.file_exists_and_not_empty(
                chimeric_genes_list_fn_full_path) == chimeric_genes_list_fn_full_path:

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


def recover_paralogs_from_all_samples_mp(list_of_sample_names,
                                         compressed_samples_set,
                                         target_genes,
                                         sampledir_parent,
                                         cpu,
                                         outdir_all_paralogs=None,
                                         outdir_no_chimeras=None):
    """

    :param list list_of_sample_names: list of sample names to process (no missing samples)
    :param compressed_samples_set: list of sample names that are *.tar.gz compressed files
    :param target_genes: list of unique gene names in the target file
    :param sampledir_parent: a path to the directory containing HybPiper output
    :param cpu: number of threads to use for multiprocessing
    :param outdir_all_paralogs: directory for output FASTA files (all paralogs)
    :param outdir_no_chimeras: directory for output FASTA files (chimeric FNA seqs skipped)
    :return:
    """

    ####################################################################################################################
    # Iterate over each sample using multiprocessing:
    ####################################################################################################################

    logger.info(f'{"[INFO]:":10} Parsing data from all samples...')

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
    # Check for chimera_check value and exit with message if not found for some samples:
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
    # Warn user if some samples didn't have a chimera check performed during `hybpiper assemble`:
    ####################################################################################################################
    samples_with_no_chimera_check_performed = [sample_name for sample_name in sample_dict_collated if not
                                               sample_dict_collated[sample_name]['chimera_check']]

    if len(samples_with_no_chimera_check_performed) != 0:
        fill = textwrap.fill(f'{"[WARNING]:":10} A chimera check was not performed during "hybpiper assemble" for '
                             f'the following samples:',
                             width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
        logger.warning(f'\n{fill}\n')

        for sample_name in sorted(list(samples_with_no_chimera_check_performed)):
            logger.warning(f'{" " * 10} {sample_name}')

        logger.warning(f'\n{" " * 10} No putative chimeric sequences will be skipped for these samples!')

    ####################################################################################################################
    # Create dictionaries to capture data for reports/heatmap:
    ####################################################################################################################
    sample_to_gene_paralog_count_dict = defaultdict(dict)
    gene_with_paralogs_to_sample_list_dict = defaultdict(list)

    ####################################################################################################################
    # Iterate over each gene and collate data from all samples:
    ####################################################################################################################
    any_putative_chimera_removed = False

    seqs_to_write_all_dict = dict()
    seqs_to_write_no_chimeras_dict = dict()

    for gene_name in target_genes:

        seqs_to_write_all_dict[gene_name] = []
        seqs_to_write_no_chimeras_dict[gene_name] = []

        ################################################################################################################
        # Iterate over each sample:
        ################################################################################################################
        for sample_name in list_of_sample_names:

            chimera_check_performed = sample_dict_collated[sample_name]['chimera_check']
            chimeric_genes_list = sample_dict_collated[sample_name]['chimeric_genes_list']
            num_seqs = 0  # default if no paralogs or *.FNA

            ############################################################################################################
            # Try to recover paralog and *.FNA seqrecords from sample dictionary:
            ############################################################################################################
            try:
                gene_paralogs_seqrecords = sample_dict_collated[sample_name]['paralog_seqrecords'][gene_name]
            except KeyError:
                gene_paralogs_seqrecords = None

            try:
                gene_fna_seqrecord = sample_dict_collated[sample_name]['fna_seqrecords'][gene_name]

                if len(gene_fna_seqrecord) != 0:
                    assert len(gene_fna_seqrecord) == 1
                    gene_fna_seqrecord = gene_fna_seqrecord[0]

            except KeyError:
                gene_fna_seqrecord = None

            if not gene_paralogs_seqrecords and not gene_fna_seqrecord:
                pass

            elif gene_paralogs_seqrecords:

                num_seqs += len(gene_paralogs_seqrecords)

                for paralog_seqrecord in gene_paralogs_seqrecords:
                    seqs_to_write_all_dict[gene_name].append(paralog_seqrecord)
                    seqs_to_write_no_chimeras_dict[gene_name].append(paralog_seqrecord)

            elif gene_fna_seqrecord:

                num_seqs = 1

                # Now skip any putative chimeric stitched contig sequences (if chimera check was performed for this
                # sample). Note that if paralogs are present they are recovered regardless of whether the gene is in
                # chimeric_genes_to_skip (which corresponds to *.FNA sequences), as paralog sequences are derived from
                # single Exonerate hits only:

                if chimera_check_performed and gene_name in chimeric_genes_list:
                    any_putative_chimera_removed = True
                    seqs_to_write_all_dict[gene_name].append(gene_fna_seqrecord)
                else:
                    seqs_to_write_no_chimeras_dict[gene_name].append(gene_fna_seqrecord)
                    seqs_to_write_all_dict[gene_name].append(gene_fna_seqrecord)

            ############################################################################################################
            # Add the number of sequences for this gene to the sample in the sample_to_gene_paralog_count_dict:
            ############################################################################################################
            sample_to_gene_paralog_count_dict[sample_name][gene_name] = num_seqs

            ############################################################################################################
            # If the gene has paralogs for this sample, add the sample name to the list for the gene in
            # gene_with_paralogs_to_sample_list_dict:
            ############################################################################################################
            if gene_paralogs_seqrecords:
                gene_with_paralogs_to_sample_list_dict[gene_name].append(sample_name)

    ####################################################################################################################
    # Set output directories:
    ####################################################################################################################
    if any_putative_chimera_removed:
        logger.info(f'{"[INFO]:":10} Creating directory: {outdir_no_chimeras}')
        if not os.path.isdir(outdir_no_chimeras):
            os.mkdir(outdir_no_chimeras)

    logger.info(f'{"[INFO]:":10} Creating directory: {outdir_all_paralogs}')
    if not os.path.isdir(outdir_all_paralogs):
        os.mkdir(outdir_all_paralogs)

    ####################################################################################################################
    # Write the 'all' and 'no chimeras' fasta files for all genes:
    ####################################################################################################################
    for gene_name, seqrecord_list in seqs_to_write_all_dict.items():
        with open(f'{outdir_all_paralogs}/{gene_name}_paralogs_all.fasta', 'w') as all_seqs_handle:
            SeqIO.write(seqrecord_list, all_seqs_handle, 'fasta')

    if any_putative_chimera_removed:  # chimera-filtered output folder will be present
        for gene_name, seqrecord_list in seqs_to_write_no_chimeras_dict.items():
            with open(f'{outdir_no_chimeras}/{gene_name}_paralogs_no_chimeras.fasta',
                      'w') as no_chimera_seqs_handle:
                SeqIO.write(seqrecord_list, no_chimera_seqs_handle, 'fasta')

    return (sample_to_gene_paralog_count_dict,
            gene_with_paralogs_to_sample_list_dict)


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

    logger.info(f'{"[INFO]:":10} Recovering paralog sequences from the HybPiper run(s)...')

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
    # Get gene names parsed from a target file.
    ####################################################################################################################
    target_genes = list(set([x.id.split('-')[-1] for x in SeqIO.parse(targetfile, 'fasta')]))

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

    ################################################################################################################
    # Check for sample names present in the namelist.txt but not found in the directory provided:
    ################################################################################################################
    list_of_sample_names = utils.check_for_missing_samples(args.namelist,
                                                           set_of_sample_names,
                                                           samples_found,
                                                           sampledir_parent,
                                                           logger)

    logger.info(f'{"[INFO]:":10} Total number of samples to process: {len(list_of_sample_names)}')

    logger.info(f'{"[INFO]:":10} Searching for paralogs for {len(list_of_sample_names)} samples, '
                f'{len(target_genes)} genes...')

    ####################################################################################################################
    # Recover sequences from all samples:
    ####################################################################################################################
    (sample_to_gene_paralog_count_dict,
     gene_with_paralogs_to_sample_list_dict) = \
        recover_paralogs_from_all_samples_mp(list_of_sample_names,
                                             compressed_samples_set,
                                             target_genes,
                                             sampledir_parent,
                                             cpu,
                                             outdir_all_paralogs=args.fasta_dir_all,
                                             outdir_no_chimeras=args.fasta_dir_no_chimeras)

    ####################################################################################################################
    # Write a *.tsv report, i.e. a matrix of sample names vs gene name with sequence counts:
    ####################################################################################################################
    with open(f'{args.paralog_report_filename}.tsv', 'w') as paralog_report_handle:
        genes = '\t'.join(target_genes)
        paralog_report_handle.write(f'Species\t{genes}\n')
        for sample, gene_paralog_count_dict in sample_to_gene_paralog_count_dict.items():
            gene_counts = '\t'.join([str(count) for count in gene_paralog_count_dict.values()])
            paralog_report_handle.write(f'{sample}\t{gene_counts}\n')

    ####################################################################################################################
    # Create a heatmap from the *.tsv file:
    ####################################################################################################################
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

    ####################################################################################################################
    # Write text statistics report:
    ####################################################################################################################
    write_paralogs_above_threshold_report(gene_with_paralogs_to_sample_list_dict,
                                          args.paralogs_list_threshold_percentage,
                                          list_of_sample_names,
                                          target_genes,
                                          args.paralogs_above_threshold_report_filename)

    logger.info(f'{"[INFO]:":10} Done!')


if __name__ == "__main__":
    standalone()
