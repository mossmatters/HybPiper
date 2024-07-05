#!/usr/bin/env python

"""
Run the assembler SPAdes with re-dos if any of the k-mers are unsuccessful.
The re-runs are attempted by removing the largest k-mer and re-running spades.
"""

import argparse
import copy
import os
import sys
import shutil
import subprocess
import re
import logging
import textwrap
from hybpiper import utils


# Create logger:
logger = logging.getLogger(f'hybpiper.assemble.{__name__}')


def make_spades_cmd_file(genelist, cov_cutoff=8, paired=True, kvals=None, unpaired=False, merged=False,
                         single_cell_mode=False):
    """
    Generates a file with gene-specific commands for running SPAdes via GNU parallel. The commands differ depending
    on options (--unpaired, --merged) and read files present for each gene.

    :param genelist: path to file containing the name of each gene that has reads distributed to its directory
    :param int cov_cutoff: coverage cutoff for SPAdes assembler
    :param bool paired: True if len(readfiles) == 2
    :param list kvals: values of k for SPAdes assemblies
    :param bool unpaired: True is an unpaired readfile has been provided for the sample
    :param bool merged: True if parameter --merged is used
    :param bool single_cell_mode: if True, run SPAdes assemblies in MDA (single-cell) mode
    :return str spades_initial_commands_file: file name of text file with SPAdes intial commands for each gene
    """

    if single_cell_mode:
        logger.info(f'{"[NOTE]:":10} Running SPAdes in MDA (single-cell) mode - be sure to check your sequences!')
        single_cell_mode_string = '--sc'
    else:
        single_cell_mode_string = ''

    if kvals:
        kvals = ','.join(kvals)

    spades_cmd_list = [f'spades.py --memory 1024 --only-assembler {single_cell_mode_string} --threads 1 --cov-cutoff',
                       str(cov_cutoff)]

    if kvals:
        spades_cmd_list.append(f'-k {kvals}')

    # Create a file of gene-specific spades commands based on options provided and read files present:
    with open(genelist, 'r') as genes_for_spades_handle:
        genes_for_spades = [gene.strip() for gene in genes_for_spades_handle]

    spades_initial_commands_file = 'spades_initial_commands.txt'
    with open(spades_initial_commands_file, 'w') as initial_spades_handle:
        for gene in genes_for_spades:
            base_list = copy.deepcopy(spades_cmd_list)
            if merged:  # implies paired read files present
                if utils.file_exists_and_not_empty(f'{gene}/{gene}_merged.fastq'):
                    base_list.append(f'--merged {gene}/{gene}_merged.fastq')
                if utils.file_exists_and_not_empty(f'{gene}/{gene}_unmerged.fastq'):
                    base_list.append(f'--12 {gene}/{gene}_unmerged.fastq')
            elif paired:
                if utils.file_exists_and_not_empty(f'{gene}/{gene}_interleaved.fasta'):
                    base_list.append(f'--12 {gene}/{gene}_interleaved.fasta')

            if utils.file_exists_and_not_empty(f'{gene}/{gene}_unpaired.fasta'):  # covers both single end or --unpaired
                base_list.append(f'-s {gene}/{gene}_unpaired.fasta')

            base_list.append(f'-o {gene}/{gene}_spades')

            full_command = ' '.join(base_list)
            initial_spades_handle.write(f'{full_command}\n')

    return spades_initial_commands_file


def spades_initial(genelist, cov_cutoff=8, cpu=None, paired=True, kvals=None, timeout=None, unpaired=False,
                   merged=False, single_cell_mode=False):
    """
    Run SPAdes on each gene separately using GNU parallel. Returns a list of genes for which the SPAdes assemblies
    failed.

    :param str genelist: path to file containing the name of each gene that has reads distributed to its directory
    :param int cov_cutoff: coverage cutoff for SPAdes assembler
    :param int cpu: number of threads/cpus to use for GNU Parallel
    :param bool paired: True if len(readfiles) == 2
    :param list kvals: values of k for SPAdes assemblies
    :param int timeout: value for GNU parallel --timeout percentage
    :param bool unpaired: True is an unpaired readfile has been provided for the sample
    :param bool merged: True if parameter --merged is used
    :param bool single_cell_mode: if True, run SPAdes assemblies in MDA (single-cell) mode
    :return: list spades_failed: list of genes for which the SPAdes assemblies failed.
    """

    if os.path.isfile("spades.log"):
        os.remove("spades.log")

    genes = [x.rstrip() for x in open(genelist)]
    logger.info(f'{"[NOTE]:":10} Running initial SPAdes assemblies for {len(genes)} genes with reads...')

    spades_initial_commands_file = make_spades_cmd_file(genelist,
                                                        cov_cutoff,
                                                        paired=paired,
                                                        kvals=kvals,
                                                        unpaired=unpaired,
                                                        merged=merged,
                                                        single_cell_mode=single_cell_mode)

    logger.info(f'{"[INFO]:":10} See file "{spades_initial_commands_file}" for a list of SPAdes commands')

    parallel_cmd_list = ['parallel', f'-j {cpu}', '--joblog', 'gnu_parallel_log.txt', '--eta']

    if timeout:
        parallel_cmd_list.append(f'--timeout {timeout}%')

    parallel_cmd_list.append(f':::: {spades_initial_commands_file} >> spades.log')

    parallel_cmd = ' '.join(parallel_cmd_list)

    try:
        result = subprocess.run(parallel_cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True, check=True)
        logger.debug(f'spades_cmd check_returncode() is: {result.check_returncode()}')
        logger.debug(f'spades_cmd stdout is: {result.stdout}')

    except subprocess.CalledProcessError as exc:
        logger.debug(f'spades_cmd FAILED. Output is: {exc}')
        logger.debug(f'spades_cmd stdout is: {exc.stdout}')
        logger.info(f'{"[WARNING]:":10} One or more genes had an error with SPAdes assembly. This may be due to '
                    f'low coverage.')

    spades_successful = []
    spades_failed = []

    for gene in genes:
        gene_failed = False
        if os.path.isfile(f'{gene}/{gene}_spades/contigs.fasta'):
            contig_file_size = os.stat(f'{gene}/{gene}_spades/contigs.fasta').st_size
            if contig_file_size > 0:
                shutil.copy(f'{gene}/{gene}_spades/contigs.fasta', f'{gene}/{gene}_contigs.fasta')
                spades_successful.append(gene)
            else:
                gene_failed = True
        else:
            gene_failed = True

        if gene_failed:
            spades_failed.append(gene)
    logger.info(f'{"[WARNING]:":10} Total number of genes with failed initial SPAdes run: {len(spades_failed)}. Gene '
                f'names can be found in the sample log file.')
    logger.debug(f'{" ".join(spades_failed)}\n')  # Write a list of genes with failed SPAdes initial runs
    return spades_failed


def rerun_spades(genelist, cov_cutoff=8, cpu=None):
    """
    Re-run SPAdes assemblies for genes that failed in the first round, removing the largest Kmer size for the failed
    run.

    :param str genelist: path to file with list of genes with failed SPAdes assemblies
    :param int cov_cutoff: coverage cutoff for SPAdes assembler
    :param int cpu: number of threads/cpus to use for GNU Parallel
    :return: list spades_duds: spades_duds is a list of genes with failed SPAdes redo assemblies.
    """

    genes = [x.rstrip() for x in open(genelist)]

    logger.info(f'{"[NOTE]:":10} Re-running SPAdes assemblies for {len(genes)} genes with unsuccessful initial '
                f'assemblies...')

    redo_cmds_file = open('redo_spades_commands.txt', 'w')

    spades_successful = []
    spades_duds = []
    genes_redos = []

    for gene in genes:
        all_kmers = [int(x[1:]) for x in os.listdir(os.path.join(gene, f'{gene}_spades')) if
                     x.startswith('K')]
        all_kmers.sort()

        if len(all_kmers) < 2:
            spades_duds.append(gene)
            continue  # i.e. don't redo the SPAdes assembly for this gene
        else:
            genes_redos.append(gene)
        redo_kmers = [str(x) for x in all_kmers[:-1]]
        restart_k = f'k{redo_kmers[-1]}'
        kvals = ','.join(redo_kmers)
        spades_cmd = f'spades.py --restart-from {restart_k} -k {kvals} --cov-cutoff {cov_cutoff} -o {gene}' \
                     f'/{gene}_spades'
        redo_cmds_file.write(spades_cmd + "\n")

    logger.info(f'{"[WARNING]:":10} In initial assemblies all Kmers failed for {len(spades_duds)} genes; these will '
                f'not be re-run')

    redo_cmds_file.close()
    redo_spades_cmd = f'parallel -j {cpu} --eta --timeout 400% :::: redo_spades_commands.txt > spades_redo.log'

    fill = textwrap.fill(f'{"[CMD]:":10} {redo_spades_cmd}', width=90, subsequent_indent=' ' * 11,
                         break_long_words=False, break_on_hyphens=False)
    logger.info(f'{fill}')

    try:
        result = subprocess.run(redo_spades_cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True,
                                check=True)
        logger.debug(f'redo_spades_cmd check_returncode() is: {result.check_returncode()}')
        logger.debug(f'redo_spades_cmd stdout is: {result.stdout}')

    except subprocess.CalledProcessError as exc:
        logger.debug(f'redo_spades_cmd FAILED. Output is: {exc}')
        logger.debug(f'redo_spades_cmd stdout is: {exc.stdout}')
        logger.info(f'{"[WARNING]:":10} One or more genes had an error with SPAdes assembly. This may be due to low '
                    f'coverage.')

    # Check whether SPAdes re-runs were successful:
    for gene in genes_redos:
        gene_failed = False
        if os.path.isfile(f'{gene}/{gene}_spades/contigs.fasta'):
            if os.stat(f'{gene}/{gene}_spades/contigs.fasta').st_size > 0:
                shutil.copy(f'{gene}/{gene}_spades/contigs.fasta', f'{gene}/{gene}_contigs.fasta')
                spades_successful.append(gene)
            else:
                gene_failed = True
        else:
            gene_failed = True

        if gene_failed:
            spades_duds.append(gene)
    logger.info(f'{"[WARNING]:":10} Total number of genes with failed SPAdes re-runs: {len(spades_duds)}. Gene names '
                f'can be found in the sample log file.')
    logger.debug(f'{" ".join(spades_duds)}\n')  # Write a list of genes with failed SPAdes re-runs
    with open('spades_duds.txt', 'w') as spades_duds_file:
        spades_duds_file.write('\n'.join(spades_duds))

    return spades_duds


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('genelist',
                        help='Text file containing the name of each gene to conduct SPAdes assembly. One gene per '
                             'line, should correspond to directories within the current directory.')
    parser.add_argument('--cpu', type=int, default=0,
                        help='Limit the number of CPUs. Default is to use all cores available.')
    parser.add_argument('--cov_cutoff', type=int, default=8, help='Coverage cutoff for SPAdes. default: %(default)s')
    parser.add_argument('--kvals', nargs='+',
                        help='Values of k for SPAdes assemblies. Default is to use SPAdes auto detection based on '
                             'read lengths (recommended).',
                        default=None)
    parser.add_argument("--redos_only", action='store_true', default=False,
                        help='Continue from previously assembled SPAdes assemblies and only conduct redos from '
                             'failed_spades.txt')
    parser.add_argument('--single', help='Reads are single end. Default is paired end.', action='store_true',
                        default=False)
    parser.add_argument("--timeout",
                        help='Use GNU Parallel to kill processes that take longer than X times the average.', default=0)
    parser.add_argument('--unpaired', help='For assembly with both paired (interleaved) and unpaired reads',
                        action="store_true", default=False)
    parser.add_argument('--merged', help='For assembly with both merged and unmerged (interleaved) reads',
                        action='store_true', default=False)
    args = parser.parse_args()

    logger.debug(f'args.merged is: {args.merged}')
    logger.debug(f'args.unpaired is: {args.unpaired}')

    if args.single:
        is_paired = False
    else:
        is_paired = True

    # Get number of cpus/threads for pipeline:
    if args.cpu:
        cpu = args.cpu
        logger.info(f'{"[INFO]:":10} Using {cpu} cpus/threads.')
    else:
        import multiprocessing
        cpu = multiprocessing.cpu_count()  # i.e. use all cpus.
        logger.info(f'{"[INFO]:":10} Number of cpus/threads not specified, using all available ({cpu}).')

    if os.path.isfile('failed_spades.txt') and args.redos_only:  # Only gets used (optional) when running standalone
        # script
        spades_failed = rerun_spades('failed_spades.txt', cpu=cpu)
    else:
        if args.unpaired:  # Create empty unpaired file if it doesn't exist
            for gene in open(args.genelist):
                gene = gene.rstrip()
                if os.path.isfile(f'{gene}/{gene}_interleaved.fasta'):
                    if not os.path.isfile(f'{gene}/{gene}_unpaired.fasta'):
                        open(f'{gene}/{gene}_unpaired.fasta', 'a').close()

        spades_failed = spades_initial(args.genelist, cov_cutoff=args.cov_cutoff, cpu=cpu, kvals=args.kvals,
                                       paired=is_paired, timeout=args.timeout, unpaired=args.unpaired,
                                       merged=args.merged)

        if len(spades_failed) > 0:
            with open('failed_spades.txt', 'w') as failed_spadefile:
                failed_spadefile.write('\n'.join(spades_failed))

            spades_duds = rerun_spades('failed_spades.txt', cov_cutoff=args.cov_cutoff, cpu=cpu)
            if len(spades_duds) == 0:
                sys.stderr.write('All redos completed successfully!\n')
            else:
                sys.exit(1)


if __name__ == '__main__':
    main()

