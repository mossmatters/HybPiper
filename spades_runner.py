#!/usr/bin/env python

"""
Run the assembler SPAdes with re-dos if any of the k-mers are unsuccessful.
The re-runs are attempted by removing the largest k-mer and re-running spades. If a final
contigs.fasta file is generated, a 'spades.ok' file is saved.
"""

import argparse
import os
import sys
import shutil
import subprocess
import re
import logging


# Create logger:
logger = logging.getLogger(f'__main__.{__name__}')


def file_exists_and_not_empty(file_name):
    """
    Check if file exists and is not empty by confirming that its size is not 0 bytes. Returns a boolean.

    :param str file_name: path to filename to check
    :return: bool
    """

    return os.path.isfile(file_name) and not os.path.getsize(file_name) == 0


def make_spades_cmd(genelist, cov_cutoff=8, cpu=None, paired=True, kvals=None, timeout=None,
                    unpaired=False, merged=False):
    """
    Generates a command string for running SPAdes via GNU parallel. Returns either a single string, or if merged=True
    returns two command strings corresponding to genes with merged and without merged sequences

    :param genelist:
    :param int cov_cutoff: coverage cutoff for SPAdes assembler
    :param int cpu: number of threads/cpus to use for GNU Parallel
    :param bool paired: True if len(readfiles) == 2
    :param list kvals: values of k for SPAdes assemblies
    :param int timeout: value for GNU parallel --timeout percentage
    :param bool unpaired: True is an unpaired readfile has been provided for the sample
    :param bool merged: True if parameter --merged is used
    :return: str spades_cmd_with_merged, spades_cmd_without_merged, OR str spades_cmd
    """

    if kvals:
        kvals = ','.join(kvals)

    parallel_cmd_list = ['time', 'parallel', '--eta']
    if cpu:
        parallel_cmd_list.append(f'-j {cpu}')
    if timeout:
        parallel_cmd_list.append(f'--timeout {timeout}%')

    spades_cmd_list = ['spades.py --only-assembler --threads 1 --cov-cutoff', str(cov_cutoff)]
    # spades_cmd_list = ["spades.py --only-assembler --sc --threads 1 --cov-cutoff", str(cov_cutoff)]  # CJJ added --sc
    if kvals:
        spades_cmd_list.append(f'-k {kvals}')
    if unpaired:
        spades_cmd_list.append('-s {}/{}_unpaired.fasta')
    if paired and not merged:
        spades_cmd_list.append('--12 {}/{}_interleaved.fasta')
    elif not merged:  # i.e. a single file of single-end reads was provided
        spades_cmd_list.append('-s {}/{}_unpaired.fasta')
    if merged:
        spades_cmd_list.append('--merged {}/{}_merged.fastq')
        spades_cmd_list.append('--12 {}/{}_unmerged.fastq')

    # Write separate gene name files for those that have merged reads and those that don't. Format the SPAdes
    # command accordingly and run as two separate subprocess commands in spades_initial():
    if merged:
        with open(genelist, 'r') as gene_file:
            contents = gene_file.readlines()
        genelist_list = [gene.strip() for gene in contents]

        genes_with_merged_reads = [gene for gene in genelist_list if file_exists_and_not_empty(
            f'{gene}/{gene}_merged.fastq')]
        with open(f'spades_genelist_with_merged.txt', 'w') as with_merged:
            for gene in genes_with_merged_reads:
                with_merged.write(f'{gene}\n')

        genes_without_merged_reads = [gene for gene in genelist_list if not file_exists_and_not_empty(
            f'{gene}/{gene}_merged.fastq')]
        with open(f'spades_genelist_without_merged.txt', 'w') as without_merged:
            for gene in genes_without_merged_reads:
                without_merged.write(f'{gene}\n')

        logger.info(f'{"[NOTE]:":10} Genes with merged reads: {len(genes_with_merged_reads)}')
        logger.info(f'{"[NOTE]:":10} Genes without merged reads: {len(genes_without_merged_reads)}')

    if merged:
        spades_cmd_list_with_merged = spades_cmd_list.copy()
        spades_cmd_list_with_merged.append('-o {}/{}_spades :::: spades_genelist_with_merged.txt >> spades.log')
        spades_cmd_list_without_merged = spades_cmd_list.copy()
        spades_cmd_list_without_merged = [re.sub('--merged {}/{}_merged.fastq', '', item) for item in
                                          spades_cmd_list_without_merged]  # Hacky - refactor
        spades_cmd_list_without_merged.append(
            '-o {}/{}_spades :::: spades_genelist_without_merged.txt >> spades.log')

        spades_cmd_with_merged = f'{" ".join(parallel_cmd_list)} {" ".join(spades_cmd_list_with_merged)}'
        spades_cmd_without_merged = f'{" ".join(parallel_cmd_list)} {" ".join(spades_cmd_list_without_merged)}'
        return spades_cmd_with_merged, spades_cmd_without_merged
    else:
        spades_cmd_list.append(f'-o {{}}/{{}}_spades :::: {genelist} > spades.log')
        spades_cmd = f'{" ".join(parallel_cmd_list)} {" ".join(spades_cmd_list)}'
        return spades_cmd


def spades_initial(genelist, cov_cutoff=8, cpu=None, paired=True, kvals=None, timeout=None, unpaired=False,
                   merged=False):
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
    :return: list spades_failed: list of genes for which the SPAdes assemblies failed.
    """

    if os.path.isfile("spades.log"):
        os.remove("spades.log")

    genes = [x.rstrip() for x in open(genelist)]

    if merged:
        spades_cmd_with_merged, spades_cmd_without_merged = make_spades_cmd(
            genelist, cov_cutoff, cpu, paired=paired, kvals=kvals, unpaired=unpaired, merged=merged, timeout=timeout)
        logger.info(f'{"[CMD:]":10} {spades_cmd_with_merged}')

        try:
            result = subprocess.run(spades_cmd_with_merged, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True)
            logger.debug(f'spades_cmd with merged check_returncode() is: {result.check_returncode()}')
            logger.debug(f'spades_cmd with merged stdout is: {result.stdout}')
            logger.debug(f'spades_cmd with merged stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            logger.debug(f'spades_cmd with merged FAILED. Output is: {exc}')
            logger.debug(f'spades_cmd with merged stdout is: {exc.stdout}')
            logger.debug(f'spades_cmd with merged stderr is: {exc.stderr}')
            logger.info(f'{"[WARN!]:":10} One or more genes had an error with SPAdes assembly. This may be due to low '
                        f'coverage. No contigs found for the following genes:\n')

        logger.info(f'[CMD]: {spades_cmd_without_merged}')
        try:
            result = subprocess.run(spades_cmd_without_merged, shell=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, universal_newlines=True)
            logger.debug(f'spades_cmd without merged check_returncode() is: {result.check_returncode()}')
            logger.debug(f'spades_cmd without merged stdout is: {result.stdout}')
            logger.debug(f'spades_cmd without merged stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            logger.debug(f'spades_cmd without merged FAILED. Output is: {exc}')
            logger.debug(f'spades_cmd without merged stdout is: {exc.stdout}')
            logger.debug(f'spades_cmd without merged stderr is: {exc.stderr}')
            logger.info(f'{"[WARN!]:":10} One or more genes had an error with SPAdes assembly. This may be due to low '
                        f'coverage. No contigs found for the following genes:\n')

    else:
        spades_cmd = make_spades_cmd(genelist, cov_cutoff, cpu, paired=paired, kvals=kvals, unpaired=unpaired,
                                     merged=merged, timeout=timeout)

        logger.info(f'{"[NOTE]:":10} Running SPAdes on {len(genes)} genes')
        logger.info(f'{"[CMD]:":10} {spades_cmd}')

        try:
            result = subprocess.run(spades_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True)
            logger.debug(f'spades_cmd check_returncode() is: {result.check_returncode()}')
            logger.debug(f'spades_cmd stdout is: {result.stdout}')
            logger.debug(f'spades_cmd stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            logger.debug(f'spades_cmd FAILED. Output is: {exc}')
            logger.debug(f'spades_cmd stdout is: {exc.stdout}')
            logger.debug(f'spades_cmd stderr is: {exc.stderr}')
            logger.info(f'{"[WARN!]:":10} One or more genes had an error with SPAdes assembly. This may be due to '
                        f'low coverage. No contigs found for the following genes:\n')

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
            logger.info(f'{" " * 11} {gene}')
            spades_failed.append(gene)
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

    redo_cmds_file = open('redo_spades_commands.txt', 'w')

    spades_successful = []
    spades_duds = []
    genes_redos = []

    for gene in genes:
        all_kmers = [int(x[1:]) for x in os.listdir(os.path.join(gene, f'{gene}_spades')) if
                     x.startswith('K')]
        all_kmers.sort()

        if len(all_kmers) < 2:
            logger.info(f'{"[WARN!]:":10} All Kmers failed for {gene}!')
            spades_duds.append(gene)
            continue
        else:
            genes_redos.append(gene)
        redo_kmers = [str(x) for x in all_kmers[:-1]]
        restart_k = f'k{redo_kmers[-1]}'
        kvals = ','.join(redo_kmers)
        spades_cmd = f'spades.py --restart-from {restart_k} -k {kvals} --cov-cutoff {cov_cutoff} -o {gene}' \
                     f'/{gene}_spades'
        redo_cmds_file.write(spades_cmd + "\n")

    redo_cmds_file.close()
    if cpu:
        redo_spades_cmd = f'parallel -j {cpu} --eta --timeout 400% :::: redo_spades_commands.txt > spades_redo.log'
    else:
        redo_spades_cmd = 'parallel --eta --timeout 400% :::: redo_spades_commands.txt > spades_redo.log'

    logger.info(f'{"[NOTE]:":10} Re-running SPAdes for {len(genes_redos)} genes\n')
    logger.info(f'{"[CMD]:":10} {redo_spades_cmd}')

    try:
        result = subprocess.run(redo_spades_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True)
        logger.debug(f'redo_spades_cmd check_returncode() is: {result.check_returncode()}')
        logger.debug(f'redo_spades_cmd stdout is: {result.stdout}')
        logger.debug(f'redo_spades_cmd stderr is: {result.stderr}')

    except subprocess.CalledProcessError as exc:
        logger.error(f'redo_spades_cmd FAILED. Output is: {exc}')
        logger.debug(f'redo_spades_cmd stdout is: {exc.stdout}')
        logger.debug(f'redo_spades_cmd stderr is: {exc.stderr}')
        logger.info(f'{"[WARN!]:":10} One or more genes had an error with SPAdes assembly. This may be due to low '
                    f'coverage. No'
                    f'contigs found for the following genes:\n')

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
            sys.stderr.write(f'{gene}\n')
            spades_duds.append(gene)
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

    if os.path.isfile('failed_spades.txt') and args.redos_only:  # Only gets used (optional) when running standalone
        # script
        spades_failed = rerun_spades('failed_spades.txt', cpu=args.cpu)
    else:
        if args.unpaired:  # Create empty unpaired file if it doesn't exist
            for gene in open(args.genelist):
                gene = gene.rstrip()
                if os.path.isfile(f'{gene}/{gene}_interleaved.fasta'):
                    if not os.path.isfile(f'{gene}/{gene}_unpaired.fasta'):
                        open(f'{gene}/{gene}_unpaired.fasta', 'a').close()

        spades_failed = spades_initial(args.genelist, cov_cutoff=args.cov_cutoff, cpu=args.cpu, kvals=args.kvals,
                                       paired=is_paired, timeout=args.timeout, unpaired=args.unpaired,
                                       merged=args.merged)

        if len(spades_failed) > 0:
            with open('failed_spades.txt', 'w') as failed_spadefile:
                failed_spadefile.write('\n'.join(spades_failed))

            spades_duds = rerun_spades('failed_spades.txt', cov_cutoff=args.cov_cutoff, cpu=args.cpu)
            if len(spades_duds) == 0:
                sys.stderr.write('All redos completed successfully!\n')
            else:
                sys.exit(1)


if __name__ == '__main__':
    main()

