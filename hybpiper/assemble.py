#!/usr/bin/env python

import argparse
import os
import sys
import shutil
import subprocess
import glob
import logging
import logging.handlers
from collections import defaultdict
import re
import textwrap
import multiprocessing
from multiprocessing import Manager
from concurrent.futures import wait, as_completed, TimeoutError, CancelledError
import time
import signal
import platform
import traceback
import pebble
import psutil

# Import HybPiper modules:
from hybpiper.version import __version__
from hybpiper import distribute_reads_to_targets
from hybpiper import distribute_targets
from hybpiper import spades_runner
from hybpiper import exonerate_hits
from hybpiper import utils

# Create logger:
logger = logging.getLogger(f'hybpiper.hybpiper_main.{__name__}')


def bwa(readfiles, targetfile, sample_dir, cpu, unpaired=False, logger=None):
    """
    Conduct a BWA search of input reads against the targetfile.

    :param str/list readfiles: list one or more read files used as input to the pipeline, or path to unpaired read file
    :param str targetfile: path to targetfile (i.e. the target file)
    :param str sample_dir: directory name for sample
    :param int cpu: number of threads/cpus to use for BWA mapping
    :param bool unpaired: True if an unpaired file has been provided, False if not
    :param logging.Logger logger: a logger object
    :return: None, or the path to the *.bam output file from BWA alignment of sample reads to the target file
    """

    targetfile_basename = os.path.basename(targetfile)

    if os.path.isfile(targetfile):
        if os.path.isfile(f'{targetfile_basename}.amb'):
            db_file = targetfile_basename
            logger.debug(f'Using existing BWA database. db_file is: {db_file}')
        else:
            logger.info(f'{"[INFO]:":10} Making nucleotide bwa index in current directory.')
            targetfiledir = os.path.split(targetfile)[0]
            if targetfiledir:
                if os.path.realpath(targetfiledir) != os.path.realpath('.'):
                    shutil.copy(targetfile, '.')
            db_file = os.path.split(targetfile)[1]
            make_bwa_index_cmd = f'bwa index {db_file}'
            fill = textwrap.fill(f'{"[CMD]:":10} {make_bwa_index_cmd}',
                                 width=90, subsequent_indent=' ' * 11, break_long_words=False, break_on_hyphens=False)
            logger.info(f'{fill}')

            try:
                result = subprocess.run(make_bwa_index_cmd,
                                        shell=True,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        universal_newlines=True,
                                        check=True)

                logger.debug(f'BWA index check_returncode() is: {result.check_returncode()}')
                logger.debug(f'BWA index stdout is: {result.stdout}')
                logger.debug(f'BWA index stderr is: {result.stderr}')

            except subprocess.CalledProcessError as exc:
                logger.error(f'BWA index FAILED. Output is: {exc}')
                logger.error(f'BWA index stdout is: {exc.stdout}')
                logger.error(f'BWA index stderr is: {exc.stderr}')
                return None

    else:
        logger.error(f'ERROR: Cannot find target file at: {targetfile}')
        return None

    if isinstance(readfiles, list) and len(readfiles) < 3:
        bwa_fastq = ' '.join(readfiles)
    elif isinstance(readfiles, str):
        bwa_fastq = readfiles  # i.e. a path (str) to a file of unpaired (different to single-end) reads
    else:
        raise ValueError(f'Can not determine whether {readfiles} is single-end, paired-end or unpaired!')

    bwa_commands = ['bwa mem', '-t', str(cpu), db_file, bwa_fastq, ' | samtools view -h -b -S - > ']
    if unpaired:
        bwa_commands.append(f'{sample_dir}_unpaired.bam')
    else:
        bwa_commands.append(f'{sample_dir}.bam')
    full_command = ' '.join(bwa_commands)
    fill = utils.fill_forward_slash(f'{"[CMD]:":10} {full_command}',
                                    width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                    break_on_forward_slash=True)
    logger.info(f'{fill}')

    try:
        result = subprocess.run(full_command,
                                shell=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                universal_newlines=True,
                                check=True)

        logger.debug(f'BWA mapping check_returncode() is: {result.check_returncode()}')
        logger.debug(f'BWA mapping stdout is: {result.stdout}')
        logger.debug(f'BWA mapping stderr is: {result.stderr}')

    except subprocess.CalledProcessError as exc:
        logger.error(f'BWA mapping FAILED. Output is: {exc}')
        logger.error(f'BWA mapping stdout is: {exc.stdout}')
        logger.error(f'BWA mapping stderr is: {exc.stderr}')
        if unpaired:
            os.remove(f'{sample_dir}_unpaired.bam')
        else:
            os.remove(f'{sample_dir}.bam')
        return None

    return f'{sample_dir}.bam'  # No return for {basename}_unpaired.bam?


def blastx(readfiles, targetfile, evalue, basename, cpu=None, max_target_seqs=10, unpaired=False, logger=None,
           diamond=False, diamond_sensitivity=False):
    """
    Creates a blast database from the full protein target file, and performs BLASTx searches of sample
    nucleotide read files against the protein database.

    :param list/str readfiles: list of paired read files OR  path to unpaired readfile if unpaired is True
    :param str targetfile: path to targetfile (i.e. the target file)
    :param float evalue: evalue to use for BLASTx searches
    :param str basename: directory name for sample
    :param int cpu: number of threads/cpus to use for BLASTx searches
    :param int max_target_seqs: maximum target sequences specified for BLASTx searches
    :param str/bool unpaired: a path if an unpaired file has been provided, boolean False if not
    :param logging.Logger logger: a logger object
    :param bool diamond: if True use DIAMOND instead of BLASTX
    :param bool/str diamond_sensitivity: sensitivity to use for DIAMOND. Default is False; uses default DIAMOND
    :return: None, or path to *.blastx output file from DIAMOND/BLASTx searches of sample reads vs targetfile
    """

    if os.path.isfile('translated_target_file.fasta'):  # i.e. a nucleotide target file was provided but not --bwa
        targetfile_basename = 'translated_target_file.fasta'
    else:  # i.e. a protein target file was provided
        targetfile_basename = os.path.basename(targetfile)

    if os.path.isfile(targetfile):
        if os.path.isfile(f'{targetfile_basename}.psq'):
            db_file = targetfile_basename
            logger.debug(f'Using existing BLAST database. db_file is: {db_file}')
        elif os.path.isfile(f'{targetfile_basename}.dmnd'):
            db_file = targetfile_basename
            logger.debug(f'Using existing DIAMOND BLAST database. db_file is: {db_file}')
        else:
            logger.info(f'{"[INFO]:":10} Making protein blastdb in current directory.')
            if os.path.split(targetfile)[0]:
                try:
                    shutil.copy(targetfile, '.')
                except shutil.SameFileError:
                    logger.debug(f'File {targetfile} already in sample directory.')
            db_file = os.path.split(targetfile)[1]
            if diamond:
                logger.info(f'{"[INFO]:":10} Using DIAMOND instead of BLASTx!')
                if diamond_sensitivity:
                    logger.info(f'{"[INFO]:":10} Using DIAMOND sensitivity "{diamond_sensitivity}"')
                makeblastdb_cmd = f'diamond makedb --in {db_file} --db {db_file}'
            else:
                makeblastdb_cmd = f'makeblastdb -dbtype prot -in {db_file}'

            fill = textwrap.fill(f'{"[CMD]:":10} {makeblastdb_cmd}', width=90, subsequent_indent=' ' * 11,
                                 break_long_words=False, break_on_hyphens=False)
            logger.info(f'{fill}')

            try:
                result = subprocess.run(makeblastdb_cmd,
                                        shell=True,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        universal_newlines=True,
                                        check=True)

                logger.debug(f'makeblastdb check_returncode() is: {result.check_returncode()}')
                logger.debug(f'makeblastdb stdout is: {result.stdout}')
                logger.debug(f'makeblastdb stderr is: {result.stderr}')

            except subprocess.CalledProcessError as exc:
                logger.error(f'makeblastdb FAILED. Output is: {exc}')
                logger.error(f'makeblastdb stdout is: {exc.stdout}')
                logger.error(f'makeblastdb stderr is: {exc.stderr}')
                return None
    else:
        logger.error(f'Cannot find target file at: {targetfile}')
        return None

    # Remove previous blast results if they exist (because we will be appending to the *.blastx output file)
    if os.path.isfile(f'{basename}.blastx'):
        os.remove(f'{basename}.blastx')

    if unpaired:
        read_file = readfiles
        # Check if read file is gzipped:
        filename, file_extension = os.path.splitext(read_file)
        if file_extension == '.gz' and not diamond:
            logger.debug(f'Processing gzipped file {os.path.basename(read_file)}')
            pipe_cmd = f"gunzip -c {read_file} | awk '{{if(NR % 4 == 1 || NR % 4 == 2) {{sub(/@/, \">\"); print; }} " \
                       f"}}'"
        elif not diamond:
            pipe_cmd = f"cat {read_file} | awk '{{if(NR % 4 == 1 || NR % 4 == 2) {{sub(/@/, \">\"); print; }} }}'"

        if diamond and diamond_sensitivity:
            blastx_command = f'diamond blastx --threads {cpu} --db {db_file} --query {read_file} --evalue {evalue} ' \
                             f'--outfmt 6 --max-target-seqs {max_target_seqs} --{diamond_sensitivity}'
        elif diamond:
            blastx_command = f'diamond blastx --threads {cpu} --db {db_file} --query {read_file} --evalue {evalue} ' \
                             f'--outfmt 6 --max-target-seqs {max_target_seqs}'
        else:
            blastx_command = f'blastx -db {db_file} -query - -evalue {evalue} -outfmt 6 -max_target_seqs' \
                             f' {max_target_seqs}'

        if not diamond:
            full_command = f"{pipe_cmd} | parallel -j {cpu} -k --block 200K --recstart '>' --pipe '{blastx_command}' " \
                           f">> {basename}_unpaired.blastx"
        else:
            full_command = f"{blastx_command} >> {basename}_unpaired.blastx"

        fill = utils.fill_forward_slash(f'{"[CMD]:":10} {full_command}',
                                        width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                        break_on_forward_slash=True)
        logger.info(f'{fill}')

        try:
            result = subprocess.run(full_command,
                                    shell=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    universal_newlines=True,
                                    check=True)

            logger.debug(f'blastx unpaired check_returncode() is: {result.check_returncode()}')
            logger.debug(f'blastx unpaired stdout is: {result.stdout}')
            logger.debug(f'blastx unpaired stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            logger.error(f'blastx unpaired FAILED. Output is: {exc}')
            logger.error(f'blastx unpaired stdout is: {exc.stdout}')
            logger.error(f'blastx unpaired stderr is: {exc.stderr}')
            return None

        return f'{basename}_unpaired.blastx'

    else:
        for read_file in readfiles:
            # Check if read file is gzipped:
            filename, file_extension = os.path.splitext(read_file)
            if file_extension == '.gz' and not diamond:
                logger.debug(f'Processing gzipped file {os.path.basename(read_file)}')
                pipe_cmd = f"gunzip -c {read_file} | awk '{{if(NR % 4 == 1 || NR % 4 == 2) {{sub(/@/, \">\"); print; " \
                           f"}} }}'"
            else:
                pipe_cmd = f"cat {read_file} | awk '{{if(NR % 4 == 1 || NR % 4 == 2) {{sub(/@/, \">\"); print; }} }}'"

            if diamond and diamond_sensitivity:
                blastx_command = f'diamond blastx --threads {cpu} --db {db_file} --query {read_file} --evalue' \
                                 f' {evalue} --outfmt 6 --max-target-seqs {max_target_seqs} --{diamond_sensitivity}'
            elif diamond:
                blastx_command = f'diamond blastx --threads {cpu} --db {db_file} --query {read_file} --evalue' \
                                 f' {evalue} --outfmt 6 --max-target-seqs {max_target_seqs}'
            else:
                blastx_command = f'blastx -db {db_file} -query - -evalue {evalue} -outfmt 6 -max_target_seqs' \
                                 f' {max_target_seqs}'

            if not diamond:
                full_command = f"{pipe_cmd} | parallel -j {cpu} -k --block 200K --recstart '>' --pipe " \
                               f"'{blastx_command}' >> {basename}.blastx"
            else:
                full_command = f"{blastx_command} >> {basename}.blastx"

            fill = utils.fill_forward_slash(f'{"[CMD]:":10} {full_command}',
                                            width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                            break_on_forward_slash=True)
            logger.info(f'{fill}')

            try:
                result = subprocess.run(full_command,
                                        shell=True,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        universal_newlines=True,
                                        check=True)

                logger.debug(f'blastx paired check_returncode() is: {result.check_returncode()}')
                logger.debug(f'blastx paired stdout is: {result.stdout}')
                logger.debug(f'blastx paired stderr is: {result.stderr}')

            except subprocess.CalledProcessError as exc:
                logger.error(f'blastx paired FAILED. Output is: {exc}')
                logger.error(f'blastx paired stdout is: {exc.stdout}')
                logger.error(f'blastx paired stderr is: {exc.stderr}')
                return None

    return f'{basename}.blastx'


def distribute_blastx(blastx_outputfile, readfiles, targetfile, target=None, unpaired_readfile=None, exclude=None,
                      merged=False, low_mem=False, logger=None):
    """
    When using blastx, distribute sample reads to their corresponding target file hits.

    :param str blastx_outputfile: tabular format output file of BLASTx search
    :param list readfiles: one or more read files used as input to the pipeline
    :param str targetfile: path to targetfile (i.e. the target file)
    :param str target: specific target(s) to use. Tab-delimited file (one <gene>\t<taxon_name> per line) or single
    taxon name
    :param str/bool unpaired_readfile: a path if an unpaired file has been provided, False if not
    :param str exclude: specify sequence not to be used as a target sequence for Exonerate
    :param bool merged: if True, write and distribute fastq files for merging with BBmerge.sh (in addition to fasta)
    :param bool low_mem: If False, reads to distribute will be saved in a dictionary and written once; uses more RAM
    :param logging.Logger logger: a logger object
    :return: None
    """

    # Distribute reads to gene directories:
    read_hit_dict = distribute_reads_to_targets.read_sorting_blastx(blastx_outputfile)

    if len(readfiles) == 2:
        logger.info(f'{"[INFO]:":10} In total, {len(read_hit_dict) * 2} reads from the paired-end read files '
                    f'will be distributed to gene directories')
        single_end = False
    elif len(readfiles) == 1:
        logger.info(f'{"[INFO]:":10} In total, {len(read_hit_dict)} reads from the single-end read file will '
                    f'be distributed to gene directories')
        single_end = True
    else:
        raise ValueError(f'Can not determine whether single-end or pair-end reads were provided!')

    distribute_reads_to_targets.distribute_reads(readfiles,
                                                 read_hit_dict,
                                                 merged=merged,
                                                 single_end=single_end,
                                                 low_mem=low_mem)

    if unpaired_readfile:
        up_blastx_outputfile = blastx_outputfile.replace('.blastx', '_unpaired.blastx')
        read_hit_dict_unpaired = distribute_reads_to_targets.read_sorting_blastx(up_blastx_outputfile)
        logger.info(f'{"[INFO]:":10} In total, {len(read_hit_dict_unpaired)} reads from the unpaired read file will be '
                    f'distributed to gene directories')

        distribute_reads_to_targets.distribute_reads([unpaired_readfile],
                                                     read_hit_dict_unpaired,
                                                     unpaired_readfile=unpaired_readfile,
                                                     low_mem=low_mem)

    # Distribute the 'best' target file sequence (translated if necessary) to each gene directory:
    if target:
        target_string = f'{target}'
    else:
        target_string = None
    if unpaired_readfile:
        unpaired_bool = True
    else:
        unpaired_bool = False
    if exclude:
        exclude_string = f'--exclude {exclude}'
    else:
        exclude_string = None

    besthits = distribute_targets.tailored_target_blast(blastx_outputfile,
                                                        unpaired_bool,
                                                        exclude_string)

    distribute_targets.distribute_targets(targetfile,
                                          delim='-',
                                          besthits=besthits,
                                          translate=False,
                                          target=target_string)
    return None


def distribute_bwa(bamfile, readfiles, targetfile, target=None, unpaired_readfile=None, exclude=None, merged=False,
                   low_mem=False, logger=None):
    """
    When using BWA mapping, distribute sample reads to their corresponding target file gene matches.

    Distribute the 'best' target file sequence (translated if the target file contains nucleotide sequences) to each
    gene directory.

    :param str bamfile: *.bam output file from BWA alignment of sample reads to the target file
    :param list readfiles: one or more read files used as input to the pipeline
    :param str targetfile: path to targetfile (i.e. the target file)
    :param str target: specific target(s) to use. Tab-delimited file (one <gene>\t<taxon_name> per line) or single
    taxon name
    :param str/bool unpaired_readfile: a path if an unpaired file has been provided, False if not
    :param str exclude: specify sequence not to be used as a target sequence for Exonerate
    :param bool merged: if True, write and distribute fastq files for merging with BBmerge.sh (in addition to fasta)
    :param bool low_mem: If False, reads to distribute will be saved in a dictionary and written once; uses more RAM
    :param logging.Logger logger: a logger object
    :return: None
    """

    # Distribute reads to gene directories:

    read_hit_dict = distribute_reads_to_targets.read_sorting_bwa(bamfile)

    if len(readfiles) == 2:
        logger.info(f'{"[INFO]:":10} In total, {len(read_hit_dict) * 2} reads from the paired-end read files '
                    f'will be distributed to gene directories')
        single_end = False

    elif len(readfiles) == 1:
        logger.info(f'{"[INFO":10} In total, {len(read_hit_dict)} reads from the single-end read file will '
                    f'be distributed to gene directories')
        single_end = True
    else:
        raise ValueError(f'Can not determine whether single-end or pair-end reads were provided!')

    distribute_reads_to_targets.distribute_reads(readfiles,
                                                 read_hit_dict,
                                                 merged=merged,
                                                 single_end=single_end,
                                                 low_mem=low_mem)

    if unpaired_readfile:
        up_bamfile = bamfile.replace('.bam', '_unpaired.bam')
        read_hit_dict_unpaired = distribute_reads_to_targets.read_sorting_bwa(up_bamfile)
        logger.info(f'{"[INFO]:":10} In total, {len(read_hit_dict_unpaired)} reads from the unpaired read file will be '
                    f'distributed to gene directories')

        distribute_reads_to_targets.distribute_reads([unpaired_readfile],
                                                     read_hit_dict_unpaired,
                                                     unpaired_readfile=unpaired_readfile,
                                                     low_mem=low_mem)

    # Distribute the 'best' target file sequence (translated if necessary) to each gene directory:
    if target:  # i.e. a target name or file of name is specified manually
        target_string = f'{target}'
    else:
        target_string = None
    if unpaired_readfile:
        unpaired_bool = True
    else:
        unpaired_bool = False
    if exclude:
        exclude_string = f'{exclude}'
    else:
        exclude_string = None

    besthits = distribute_targets.tailored_target_bwa(bamfile,
                                                      unpaired_bool,
                                                      exclude_string)

    distribute_targets.distribute_targets(targetfile,
                                          delim='-',
                                          besthits=besthits,
                                          translate=True,
                                          target=target_string)
    return None


def spades(genes, cov_cutoff=8, cpu=None, paired=True, kvals=None, timeout=None, unpaired=False,
           merged=False, logger=None, keep_folder=False, single_cell_mode=False):
    """
    Run SPAdes on each gene separately using GNU parallel.

    :param list genes: a list of genes names that have reads distributed to their directories
    :param int cov_cutoff: coverage cutoff for SPAdes assembler
    :param int cpu: number of threads/cpus to use for GNU Parallel
    :param bool paired: True if len(readfiles) == 2
    :param list kvals: values of k for SPAdes assemblies
    :param int timeout: value for GNU parallel --timeout percentage
    :param bool unpaired: True is an unpaired readfile has been provided for the sample, else False
    :param bool merged: True if parameter --merged is used, else False
    :param logging.Logger logger: a logger object
    :param bool keep_folder: if True, don't delete the SPAdes assembly folder after contig recovery
    :param bool single_cell_mode: if True, run SPAdes assemblies in MDA (single-cell) mode
    :return: list spades_genelist: a list of gene names that had successful SPAdes assemblies (contigs.fasta produced)
    """

    with open('spades_genelist.txt', 'w') as spadesfile:
        spadesfile.write('\n'.join(genes) + '\n')

    if os.path.isfile('spades.log'):
        os.remove('spades.log')
    if os.path.isfile('spades_redo.log'):
        os.remove('spades_redo.log')

    logger.debug(f'args.merged is: {merged}')
    logger.debug(f'args.unpaired is: {unpaired}')

    spades_failed = spades_runner.spades_initial('spades_genelist.txt',
                                                 cov_cutoff=cov_cutoff,
                                                 cpu=cpu,
                                                 kvals=kvals,
                                                 paired=paired,
                                                 timeout=timeout,
                                                 unpaired=unpaired,
                                                 merged=merged,
                                                 single_cell_mode=single_cell_mode)

    logger.info(f'{"[INFO]:":10} Finished running initial SPAdes assemblies for all genes with reads!')
    if len(spades_failed) > 0:
        with open('failed_spades.txt', 'w') as failed_spadefile:
            failed_spadefile.write('\n'.join(spades_failed))

        spades_duds = spades_runner.rerun_spades('failed_spades.txt',
                                                 cov_cutoff=cov_cutoff,
                                                 cpu=cpu)

        logger.info(f'{"[INFO]:":10} Finished re-running SPAdes assemblies for genes with unsuccessful initial '
                    f'assemblies!')

        if len(spades_duds) == 0:
            logger.info(f'{"[INFO]:":10} All SPAdes re-runs completed successfully!')

    if os.path.isfile('spades_duds.txt'):  # Written by spades_runner.rerun_spades()
        spades_duds = [x.rstrip() for x in open('spades_duds.txt')]
    else:
        spades_duds = []

    spades_genelist = []
    for gene in genes:
        if gene not in set(spades_duds):
            spades_genelist.append(gene)
        if not keep_folder:
            logger.debug(f'Deleting SPAdes assembly folder for gene {gene}.')
            if os.path.isdir(f'{gene}/{gene}_spades'):
                shutil.rmtree(f'{gene}/{gene}_spades')

    with open('exonerate_genelist.txt', 'w') as genefile:
        genefile.write('\n'.join(spades_genelist) + '\n')

    return spades_genelist


def exonerate(gene_name,
              sample_dir,
              pid_list,
              thresh=55,
              paralog_warning_min_length_percentage=0.75,
              depth_multiplier=10,
              no_stitched_contig=False,
              stitched_contig_pad_n=True,
              chimera_check=False,
              bbmap_memory=250,
              bbmap_subfilter=7,
              bbmap_threads=2,
              chimeric_stitched_contig_edit_distance=5,
              chimeric_stitched_contig_discordant_reads_cutoff=5,
              worker_configurer_func=None,
              counter=None,
              lock=None,
              genes_to_process=0,
              no_intronerate=False,
              no_padding_supercontigs=False,
              keep_intermediate_files=False,
              exonerate_hit_sliding_window_size=3,
              exonerate_hit_sliding_window_thresh=55,
              exonerate_skip_frameshifts=False,
              exonerate_skip_internal_stops=False,
              exonerate_skip_terminal_stops=False,
              verbose_logging=False):
    """
    :param str gene_name: name of a gene that had at least one SPAdes contig
    :param str sample_dir: directory name for sample
    :param multiprocessing.managers.ListProxy pid_list: list shared by processes for capturing parent PIDs
    :param int thresh: percent identity threshold for stitching together Exonerate results
    :param float paralog_warning_min_length_percentage: min % of a contig vs ref protein length for a paralog warning
    :param int depth_multiplier: assign long paralog as main if coverage depth <depth_multiplier> other paralogs
    :param bool no_stitched_contig: if True, don't create stitched contigs and just use longest Exonerate hit
    :param bool stitched_contig_pad_n: if True, pad gaps in stitched contig with Ns corresponding to query gap * 3
    :param bool chimera_check: run chimera check. Default is False
    :param int bbmap_memory: MB memory (RAM ) to use for bbmap.sh
    :param int bbmap_subfilter: ban alignments with more than this many substitutions
    :param int bbmap_threads: number of threads to use for BBmap when searching for chimeric stitched contigs
    :param int chimeric_stitched_contig_edit_distance: min num differences for a read pair to be flagged as discordant
    :param int chimeric_stitched_contig_discordant_reads_cutoff: min num discordant reads pairs to flag a stitched
    contig as chimeric
    :param function worker_configurer_func: function to configure logging to file
    :param multiprocessing.managers.ValueProxy counter:
    :param multiprocessing.managers.AcquirerProxy lock:
    :param int genes_to_process: total number of genes to be processed via Exonerate
    :param bool no_intronerate: if True, do not run intronerate
    :param bool no_padding_supercontigs: if True, don't pad contig joins in supercontigs with stretches if 10 Ns
    :param bool keep_intermediate_files: if True, keep intermediate files from stitched contig and intronerate()
    processing
    :param int exonerate_hit_sliding_window_size: size of the sliding window (in amino-acids) when trimming termini
    of Exonerate hits
    :param int exonerate_hit_sliding_window_thresh: percentage similarity threshold for the sliding window (in
    amino-acids) when trimming termini of Exonerate hits
    :param bool exonerate_skip_frameshifts: skip Exonerate hits where SPAdes sequence contains frameshifts
    :param bool exonerate_skip_internal_stops: skip Exonerate hits where SPAdes sequence contains internal stop codons
    :param bool exonerate_skip_terminal_stops: skip Exonerate hits where SPAdes sequence contains a terminal stop codon
    :param bool verbose_logging: if True, log additional information to file
    :return: str gene_name, str prot_length OR None, None
    """

    # get parent PID and add to shared list; this is used to kill child processes on user interrupt:
    pid = os.getpid()
    pid_list.append(pid)

    logger = logging.getLogger()  # Assign root logger from inside the new Python process (ProcessPoolExecutor pool)
    if logger.hasHandlers():
        logger.handlers.clear()

    worker_configurer_func(gene_name)  # set up process-specific logging to file
    logger = logging.getLogger(gene_name)
    logger.setLevel(logging.DEBUG)

    # Write gene name, start time, PID etc. to a dictionary shared by the multiprocessing pool:
    start = time.time()
    worker_stats = (gene_name, multiprocessing.current_process().pid, multiprocessing.Process().name, start)
    if verbose_logging:
        logger.debug(f'worker_stats for {gene_name} are: {worker_stats}')

    # Create directories for output files based on the prefix name, or assemblyfile name:
    prefix = exonerate_hits.create_output_directories(f'{gene_name}/{sample_dir}',
                                                      f'{gene_name}/{gene_name}_contigs.fasta')
    logger.debug(f'prefix is: {prefix}')

    # Set whether the chimeric stitched contig test will be performed, and whether a file of interleaved reads is found:
    perform_stitched_contig_chimera_test, path_to_interleaved_fasta = exonerate_hits.set_stitched_contig_chimera_test(
        no_stitched_contig,
        prefix)

    logger.debug(f'perform_stitched_contig_chimera_test is: {perform_stitched_contig_chimera_test}')
    logger.debug(f'path_to_interleaved_fasta is: {path_to_interleaved_fasta}')

    # Read the SPAdes contigs and the 'best' protein reference seq into SeqIO dictionaries:
    try:
        spades_assembly_dict, best_protein_ref_dict = exonerate_hits.parse_spades_and_best_reference(
            f'{gene_name}/{gene_name}_contigs.fasta',
            f'{gene_name}/{gene_name}_target.fasta',
            prefix)

        if verbose_logging:
            logger.debug(f'spades_assembly_dict is: {spades_assembly_dict}')
            logger.debug(f'best_protein_ref_dict is: {best_protein_ref_dict}')

    except FileNotFoundError as e:
        logger.error(f"\n{'[ERROR!]:':10} Couldn't find an expected file for either the SPAdes assembly or the protein "
                     f"reference for gene {gene_name}, error is {e}")
        with lock:
            counter.value += 1
            sys.stdout.write(f'\r{"[INFO]:":10} Finished running Exonerate for gene {gene_name}, {counter.value}'
                             f'/{genes_to_process}')

        end = time.time()
        proc_run_time = end - start
        return gene_name, None, proc_run_time  # return gene_name to that log can be re-logged to main log file

    # Set variable for intronerate_success
    intronerate_success = 'N/A'

    # Perform Exonerate search with 'best' protein ref as query and SPAdes contigs as subjects:
    exonerate_text_output = exonerate_hits.initial_exonerate(f'{gene_name}/{gene_name}_target.fasta',
                                                             f'{gene_name}/{gene_name}_contigs.fasta',
                                                             prefix)

    if exonerate_text_output:  # i.e. if the initial_exonerate DID produce a result
        exonerate_result = exonerate_hits.parse_exonerate_and_get_stitched_contig(
            exonerate_text_output,
            query_file=f'{gene_name}/{gene_name}_target.fasta',
            paralog_warning_min_length_percentage=paralog_warning_min_length_percentage,
            thresh=thresh,
            logger=logger,
            prefix=prefix,
            chimera_check=chimera_check,
            discordant_cutoff=chimeric_stitched_contig_discordant_reads_cutoff,
            edit_distance=chimeric_stitched_contig_edit_distance,
            bbmap_subfilter=bbmap_subfilter,
            bbmap_memory=bbmap_memory,
            bbmap_threads=bbmap_threads,
            interleaved_fasta_file=path_to_interleaved_fasta,
            no_stitched_contig=no_stitched_contig,
            stitched_contig_pad_n=stitched_contig_pad_n,
            spades_assembly_dict=spades_assembly_dict,
            depth_multiplier=depth_multiplier,
            keep_intermediate_files=keep_intermediate_files,
            exonerate_hit_sliding_window_size=exonerate_hit_sliding_window_size,
            exonerate_hit_sliding_window_thresh=exonerate_hit_sliding_window_thresh,
            exonerate_skip_frameshifts=exonerate_skip_frameshifts,
            exonerate_skip_internal_stops=exonerate_skip_internal_stops,
            exonerate_skip_terminal_stops=exonerate_skip_terminal_stops,
            verbose_logging=verbose_logging)

        if not no_intronerate and exonerate_result and exonerate_result.hits_filtered_by_pct_similarity_dict:
            if verbose_logging:
                logger.debug(f'exonerate_result.hits_subsumed_hits_removed_overlaps_trimmed_dict for gene {gene_name} '
                             f'is: {exonerate_result.hits_subsumed_hits_removed_overlaps_trimmed_dict}')
            logger.debug(f'Running intronerate')

            try:
                exonerate_hits.intronerate(exonerate_result,
                                           spades_assembly_dict,
                                           logger=logger,
                                           no_padding_supercontigs=no_padding_supercontigs,
                                           keep_intermediate_files=keep_intermediate_files,
                                           verbose_logging=verbose_logging)

                intronerate_success = True

            except Exception:
                tb = traceback.format_exc()
                logger.debug(f'Intronerate failed for gene {gene_name} with Exception:\n{tb}')
                intronerate_success = False
    else:
        exonerate_result = False

    with lock:
        counter.value += 1
        sys.stdout.write(f'\r{"[INFO]:":10} Finished running Exonerate for gene {gene_name}, {counter.value}'
                         f'/{genes_to_process}')

    if not exonerate_text_output or not exonerate_result or not exonerate_result.stitched_contig_seqrecord:
        end = time.time()
        proc_run_time = end - start

        # return gene_name to that exonerate_hits.py log can be re-logged to main log file:
        return gene_name, None, proc_run_time, False, intronerate_success

    end = time.time()
    proc_run_time = end - start

    return gene_name, \
           len(exonerate_result.stitched_contig_seqrecord), \
           proc_run_time, \
           exonerate_result.stop_codons_in_seqrecord_bool, \
           intronerate_success


def exonerate_multiprocessing(genes,
                              sample_dir,
                              thresh=55,
                              paralog_warning_min_length_percentage=0.75,
                              pool_threads=None,
                              depth_multiplier=10,
                              no_stitched_contig=False,
                              stitched_contig_pad_n=True,
                              chimera_check=False,
                              bbmap_memory=250,
                              bbmap_subfilter=7,
                              bbmap_threads=2,
                              chimeric_stitched_contig_edit_distance=5,
                              chimeric_stitched_contig_discordant_reads_cutoff=5,
                              logger=None,
                              no_intronerate=False,
                              no_padding_supercontigs=False,
                              keep_intermediate_files=False,
                              exonerate_contigs_timeout=None,
                              exonerate_hit_sliding_window_size=3,
                              exonerate_hit_sliding_window_thresh=55,
                              exonerate_skip_frameshifts=False,
                              exonerate_skip_internal_stops=False,
                              exonerate_skip_terminal_stops=False,
                              verbose_logging=False):
    """
    Runs the function exonerate() using multiprocessing.

    :param list genes: list of genes that had successful SPAdes runs
    :param str sample_dir: directory name for sample
    :param int thresh: percent identity threshold for stitching together Exonerate results
    :param float paralog_warning_min_length_percentage: min % of a contig vs ref protein length for a paralog warning
    :param int pool_threads: number of threads/cpus to use for the pebble.ProcessPool pool
    :param int depth_multiplier: assign long paralog as main if coverage depth <depth_multiplier> other paralogs
    :param bool no_stitched_contig: if True, don't create stitched contig and just use longest Exonerate hit
    :param bool stitched_contig_pad_n: if True, pad gaps in stitched contig with Ns corresponding to query gap * 3
    :param bool chimera_check: run chimera check. Default is False
    :param int bbmap_memory: MB memory (RAM ) to use for bbmap.sh
    :param int bbmap_subfilter: ban alignments with more than this many substitutions
    :param int bbmap_threads: number of threads to use for BBmap when searching for chimeric stitched contigs
    :param int chimeric_stitched_contig_edit_distance: min num differences for a read pair to be flagged as discordant
    :param int chimeric_stitched_contig_discordant_reads_cutoff: min num discordant reads pairs to flag a stitched
    contig as chimeric
    :param logging.Logger logger: a logger object
    :param bool no_intronerate: if True, intronerate will not be run (if a gene is constructed from hits with introns)
    :param bool no_padding_supercontigs: if True, don't pad contig joins in supercontigs with stretches if 10 Ns
    :param bool keep_intermediate_files: if True, keep individual Exonerate logs rather than deleting them after
    re-logging to the main sample log file
    :param int exonerate_contigs_timeout: number of second for pebble.ProcessPool pool.schedule timeout
    :param int exonerate_hit_sliding_window_size: size of the sliding window (in amino-acids) when trimming termini
    of Exonerate hits
    :param int exonerate_hit_sliding_window_thresh: percentage similarity threshold for the sliding window (in
    amino-acids) when trimming termini of Exonerate hits
    :param bool exonerate_skip_frameshifts: skip Exonerate hits where SPAdes sequence contains frameshifts
    :param bool exonerate_skip_internal_stops: skip Exonerate hits where SPAdes sequence contains internal stop codons
    :param bool exonerate_skip_terminal_stops: skip Exonerate hits where SPAdes sequence contains a terminal stop codon
    :param bool verbose_logging: if True, log additional information to file
    :return:
    """

    logger.debug(f'no_intronerate is: {no_intronerate}')
    logger.debug(f'exonerate_contigs_timeout is: {exonerate_contigs_timeout}')
    logger.debug(f'exonerate_hit_sliding_window_size is: {exonerate_hit_sliding_window_size}')
    logger.debug(f'exonerate_hit_sliding_window_thresh is: {exonerate_hit_sliding_window_thresh}')
    logger.debug(f'exonerate_skip_frameshifts is: {exonerate_skip_frameshifts}')
    logger.debug(f'exonerate_skip_internal_stops is: {exonerate_skip_internal_stops}')
    logger.debug(f'exonerate_skip_terminal_stops is: {exonerate_skip_terminal_stops}')
    logger.debug(f'chimera_check is: {chimera_check}')
    logger.debug(f'stitched_contig_pad_n is: {stitched_contig_pad_n}')
    logger.debug(f'exonerate_multiprocessing pool_threads is: {pool_threads}')

    logger.info(f'{"[INFO]:":10} Running exonerate_hits for {len(genes)} genes...')
    genes_to_process = len(genes)

    try:
        with pebble.ProcessPool(max_workers=pool_threads, context=multiprocessing.get_context('fork')) as pool:
            genes_cancelled_due_to_timeout = []
            genes_cancelled_due_to_errors = []
            future_results_dict = defaultdict()
            manager = Manager()
            lock = manager.Lock()
            pid_list = manager.list()
            counter = manager.Value('i', 0)
            kwargs_for_schedule = {"thresh": thresh,
                                   "paralog_warning_min_length_percentage": paralog_warning_min_length_percentage,
                                   "depth_multiplier": depth_multiplier,
                                   "no_stitched_contig": no_stitched_contig,
                                   "stitched_contig_pad_n": stitched_contig_pad_n,
                                   "chimera_check": chimera_check,
                                   "bbmap_memory": bbmap_memory,
                                   "bbmap_subfilter": bbmap_subfilter,
                                   "bbmap_threads": bbmap_threads,
                                   "chimeric_stitched_contig_edit_distance": chimeric_stitched_contig_edit_distance,
                                   "chimeric_stitched_contig_discordant_reads_cutoff":
                                       chimeric_stitched_contig_discordant_reads_cutoff,
                                   "worker_configurer_func": utils.worker_configurer,
                                   "counter": counter,
                                   "lock": lock,
                                   "genes_to_process": genes_to_process,
                                   "no_intronerate": no_intronerate,
                                   "no_padding_supercontigs": no_padding_supercontigs,
                                   "keep_intermediate_files": keep_intermediate_files,
                                   "exonerate_hit_sliding_window_size": exonerate_hit_sliding_window_size,
                                   "exonerate_hit_sliding_window_thresh": exonerate_hit_sliding_window_thresh,
                                   "exonerate_skip_frameshifts": exonerate_skip_frameshifts,
                                   "exonerate_skip_internal_stops": exonerate_skip_internal_stops,
                                   "exonerate_skip_terminal_stops": exonerate_skip_terminal_stops,
                                   "verbose_logging": verbose_logging}

            for gene_name in genes:  # schedule jobs and store each future in a future : gene_name dict
                exonerate_job = pool.schedule(exonerate, args=[gene_name, sample_dir, pid_list],
                                              kwargs=kwargs_for_schedule, timeout=exonerate_contigs_timeout)
                future_results_dict[exonerate_job] = gene_name

            futures_list = [future for future in future_results_dict.keys()]

            # As per-gene Exonerate runs complete, read the gene log, log it to the main logger, delete gene log:
            genes_with_failed_intronerate = []

            with open('genes_with_seqs.txt', 'w') as genes_with_seqs_handle:
                with open(f'{sample_dir}_genes_with_non_terminal_stop_codons.txt', 'w') as genes_with_stops_handle:
                    for future in as_completed(futures_list):
                        try:
                            gene_name, prot_length, run_time, stop_codons_in_seqrecord_bool, intronerate_success = \
                                future.result()
                            if gene_name:  # i.e. log the Exonerate run regardless of success
                                gene_log_file_list = glob.glob(f'{gene_name}/{gene_name}*log')
                                gene_log_file_list.sort(key=os.path.getmtime)  # sort by time in case of previous undeleted log
                                gene_log_file_to_cat = gene_log_file_list[-1]  # get most recent gene log
                                with open(gene_log_file_to_cat) as gene_log_handle:
                                    lines = gene_log_handle.readlines()
                                    for line in lines:
                                        logger.debug(line.strip())  # log contents to main logger
                                if not keep_intermediate_files:
                                    os.remove(gene_log_file_to_cat)  # delete the Exonerate log file

                            # Write the 'gene_name', 'prot_length' strings returned by each process to file:
                            if gene_name and prot_length:
                                genes_with_seqs_handle.write(f'{gene_name}\t{prot_length}\n')

                            if stop_codons_in_seqrecord_bool:
                                genes_with_stops_handle.write(f'{gene_name}\n')

                            # Check if intronerate was run and, if so, if it ran successfully:
                            if intronerate_success == 'N/A':
                                pass
                            elif intronerate_success:
                                pass
                            else:
                                genes_with_failed_intronerate.append(future_results_dict[future])

                        except TimeoutError as err:
                            logger.debug(f'\nProcess timeout - exonerate() for gene {future_results_dict[future]} took '
                                         f'more than {err.args[1]} seconds to complete and was cancelled')
                            genes_cancelled_due_to_timeout.append(future_results_dict[future])
                        except CancelledError:
                            logger.debug(f'CancelledError raised for gene {future_results_dict[future]}')
                        except Exception as error:
                            genes_cancelled_due_to_errors.append(future_results_dict[future])
                            logger.debug(f'For gene {future_results_dict[future]} exonerate() raised: {error}')
                            tb = traceback.format_exc()
                            logger.debug(f'traceback is:\n{tb}')

        wait(futures_list, return_when="ALL_COMPLETED")  # redundant, but...

        if genes_with_failed_intronerate:
            fill = textwrap.fill(f'{"[WARNING]:":10} The Intronerate step of the pipeline failed for the '
                                 f'following genes:\n', width=90, subsequent_indent=" " * 11)
            logger.info('')
            logger.info(fill)
            for gene in genes_with_failed_intronerate:
                logger.info(f'{" " * 11}{gene}')

            logger.info(f'\nPlease see the log file in the sample directory for more information.')

        if genes_cancelled_due_to_errors:
            fill = textwrap.fill(f'{"[WARNING]:":10} The exonerate_contigs step of the pipeline failed for the '
                                 f'following genes:\n', width=90, subsequent_indent=" " * 11)
            logger.info('')
            logger.info(fill)
            for gene in genes_cancelled_due_to_errors:
                logger.info(f'{" " * 11}{gene}')

            logger.info(f'\nPlease see the log file in the sample directory for more information.')

        if genes_cancelled_due_to_timeout:
            fill = textwrap.fill(f'{"[WARNING]:":10} The exonerate_contigs step of the pipeline was cancelled for the '
                                 f'following genes, due to exceeding the timeout limit of {exonerate_contigs_timeout} '
                                 f'seconds\n:', width=90, subsequent_indent=" " * 11)
            logger.info('')
            logger.info(fill)
            for gene in genes_cancelled_due_to_timeout:
                logger.info(f'{" " * 11}{gene}')

            fill = textwrap.fill(f'{"[INFO]:":10} This is most likely caused by many low-complexity reads mapping to '
                                 f'the corresponding gene sequences in the target file, resulting in a SPAdes assembly '
                                 f'with many (i.e. hundreds) of repetitive and low-complexity contigs. Subsequently, '
                                 f'Exonerate searches of these many low-complexity contigs can take a long time. We '
                                 f'strongly recommend removing such low-complexity sequences from your target file. '
                                 f'The command "hybpiper check_targetfile" can assist in identifying these '
                                 f'sequences.', width=90, subsequent_indent=" " * 11)
            logger.info(fill)

    except KeyboardInterrupt:
        signal.signal(signal.SIGINT, signal.SIG_IGN)  # Ignore additional SIGINT while HybPiper cleans up
        pid_set = set(pid_list)
        parent_list = []
        for process_pid in pid_set:
            parent = psutil.Process(process_pid)
            parent_list.append(parent)
        logger.info(f'\n\nExiting HybPiper due to user interrupt, please wait a moment...\n')
        while True:
            count = 0
            child_list = [parent.children(recursive=True) for parent in parent_list]
            for children_list in child_list:
                for child in children_list:
                    count += 1
                    try:
                        child.kill()
                    except psutil.NoSuchProcess:
                        pass
            if count == 0:
                break
        sys.exit(1)


def main(args):
    """
    Assemble gene, intron, supercontig and paralog sequences.

    :param argparse.Namespace args: argparse namespace with subparser options for function main()
    :return None: no return value specified; default is None
    """

    # Check that args.prefix is NOT a path:
    if args.prefix and re.search('/', args.prefix):
        sys.exit(f'{"[ERROR]:":10} The prefix supplied with the parameter --prefix must not contain forward '
                 f'slashes. You supplied: {args.prefix}')

    # Get a list of read files from args.readfiles (doesn't include any read file passed in via --unpaired flag):
    readfiles = [os.path.abspath(x) for x in args.readfiles]

    if len(readfiles) > 2:
        sys.exit(f'{"[ERROR]:":10} Please provide a maximum of two read files (R1 and R2) to the -r / '
                 f'--readfiles parameter')

    # Generate a directory for the sample:
    parent_dir, sample_dir = utils.make_basename(args.readfiles,
                                                 prefix=args.prefix,
                                                 output_folder=args.output_folder)

    full_sample_directory = os.path.join(parent_dir, sample_dir)

    # Create logger:
    logger = utils.setup_logger(__name__, f'{full_sample_directory}/{sample_dir}_hybpiper_assemble')

    # Log command line and default parameters:
    logger.info(f'{"[INFO]:":10} HybPiper version {__version__} was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv),
                         width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11, break_on_hyphens=False)
    logger.info(f'{fill}\n')
    logger.debug(args)

    # Log system details for debugging:
    utils.get_platform_info(logger=logger)

    if platform.system() == 'Darwin':
        utils.check_macos_version(logger=logger)

    # Log ulimit details for debugging:
    utils.get_ulimit_info(logger=logger)

    # Get number of cpus/threads for pipeline:
    if args.cpu:
        cpu = args.cpu
        logger.info(f'{"[INFO]:":10} Using {cpu} cpus/threads.')
    else:
        cpu = multiprocessing.cpu_count() - 1  # i.e. use all cpus.
        logger.info(f'{"[INFO]:":10} Number of cpus/threads not specified, using all available cpus minus 1 ({cpu}).')

    ####################################################################################################################
    # Check dependencies
    ####################################################################################################################
    if utils.check_dependencies(logger=logger):
        logger.info(f'{"[INFO]:":10} Everything looks good!')
    else:
        logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
        sys.exit()

    ####################################################################################################################
    # Check that --start_from <= --end_with
    ####################################################################################################################
    assemble_stages_dict = {'map_reads': 0,
                            'distribute_reads': 1,
                            'assemble_reads': 2,
                            'exonerate_contigs': 3}

    if not assemble_stages_dict[args.start_from] <= assemble_stages_dict[args.end_with]:
        logger.error(f'{"[ERROR]:":10} The selected "--start_from" step is greater than the selected "--end_with" '
                     f'step!')
        logger.error(f'{" " * 10} --start_from: {args.start_from}')
        logger.error(f'{" " * 10} --end_with: {args.end_with}')
        logger.error('')
        logger.error(f'{" " * 10} The order of steps is: map_reads, distribute_reads, assemble_reads, '
                     f'exonerate_contigs')
        sys.exit()

    ####################################################################################################################
    # Check read and target files
    ####################################################################################################################
    # Set target file type and path, and check it exists and isn't empty:
    if args.targetfile_dna:
        targetfile = os.path.abspath(args.targetfile_dna)
        targetfile_type = 'DNA'
    elif args.targetfile_aa:
        targetfile = os.path.abspath(args.targetfile_aa)
        targetfile_type = 'protein'

    if os.path.isfile(targetfile) and not os.path.getsize(targetfile) == 0:
        logger.debug(f'Input target file {os.path.basename(targetfile)} exists and is not empty, proceeding...')
    else:
        sys.exit(f'Input target file {os.path.basename(targetfile)} does not exist or is empty!')

    logger.debug(f'The target file {os.path.basename(targetfile)} has been provided, containing {targetfile_type} '
                 f'sequences')

    # Check that the input read files exist and aren't empty:
    for read_file in readfiles:
        if os.path.isfile(read_file) and not os.path.getsize(read_file) == 0:
            logger.debug(f'Input read file {read_file} exists and is not empty, proceeding...')
        else:
            sys.exit(f'Input read file {read_file} does not exist or is empty!')

    if args.unpaired:
        unpaired_readfile = os.path.abspath(args.unpaired)
        if os.path.isfile(unpaired_readfile) and not os.path.getsize(unpaired_readfile) == 0:
            logger.debug(f'Input read file {os.path.basename(unpaired_readfile)} exists and is not empty, '
                         f'proceeding...')
        else:
            sys.exit(f'Input read file {os.path.basename(unpaired_readfile)} does not exist or is empty!')
    else:
        unpaired_readfile = False

    # If only a single readfile is supplied, set --merged to False regardless of user input:
    if len(readfiles) == 1 and args.merged:
        logger.info(f'{"[INFO]:":10} The flag --merged has been provided but only a single read file has been '
                    f'supplied. Setting --merged to False.')
        args.merged = False

    # If a file of unpaired reads is provided via the --unpaired parameter and only a single readfile is provided via
    # the -r/--readfiles parameter, exit with an error message:
    if len(readfiles) == 1 and args.unpaired:
        sys.exit(f'{"[ERROR]:":10} You have provided a single file of reads using the -r/--readfiles parameter ('
                 f'{os.path.basename(readfiles[0])}), along with a file of unpaired reads via the --unpaired '
                 f'parameter ({os.path.basename(args.unpaired)}). Please concatenate these two files and provide the '
                 f'single file as input using the -r/--readfiles parameter')

    # Check that the target file is formatted correctly and translates correctly. If it contains DNA sequences but
    # arg.bwa is false, translate and return a list of translated sequences instead of file path:
    targetfile = utils.check_targetfile(targetfile,
                                        targetfile_type,
                                        full_sample_directory=full_sample_directory,
                                        using_bwa=args.bwa,
                                        skip_targetfile_checks=args.skip_targetfile_checks,
                                        logger=logger)

    ####################################################################################################################
    # Check manually provided targets if provided via the parameter --target
    ####################################################################################################################
    if args.target:
        target_path = os.path.abspath(args.target)
        if os.path.isfile(target_path):
            logger.debug(f'A file of preferred target sequences for Exonerate searches has been provided: '
                         f'{target_path}')
            target = target_path
        else:
            logger.debug(f'A single preferred target sequence taxon name for Exonerate searches has been provided: '
                         f'{args.target}')
            target = args.target
    else:
        target = None

    ####################################################################################################################
    # Check if files exist from a previous run and require --force_overwrite of so:
    ####################################################################################################################
    files_from_previous_run_dict = utils.check_for_previous_run_output(full_sample_directory,
                                                                       args.start_from,
                                                                       args.end_with,
                                                                       assemble_stages_dict)

    if len(files_from_previous_run_dict) != 0:
        fill = textwrap.fill(f'{"[WARNING]:":10} Output files from a previous run have been detected for the '
                             f'selected pipeline steps: "--start_from {args.start_from}" to "--end_with '
                             f'{args.end_with}". Selected steps with existing output files are:',
                             width=90, subsequent_indent=' ' * 11)
        logger.info(fill)
        logger.info('')
        for step in files_from_previous_run_dict.keys():
            logger.info(f'{" " * 10} {step}')
        logger.info('')

        if args.force_overwrite:
            fill = textwrap.fill(f'{"[WARNING]:":10} Option "--force_overwrite" provided - overwriting all '
                                 f'existing output files for selected pipelines steps: "{args.start_from}" to '
                                 f'"{args.end_with}"!',
                                 width=90, subsequent_indent=' ' * 11)
            logger.info(fill)
        else:
            fill = textwrap.fill(f'{"[ERROR]:":10} Please provide option "--force_overwrite" if you want to '
                                 f'overwrite all existing output files for selected pipelines steps!',
                                 width=90, subsequent_indent=' ' * 11)
            logger.info(fill)
            sys.exit()

    # Log output folder:
    logger.info(f'{"[INFO]:":10} Output will be written to the directory: {full_sample_directory}')

    # Move in to the sample directory:
    os.chdir(full_sample_directory)

    ####################################################################################################################
    # Map reads to target file sequences
    ####################################################################################################################

    if args.start_from not in ['map_reads']:  # 'map_reads' is the default
        logger.info(f'{"[INFO]:":10} Parameter "--start_from {args.start_from}" supplied, skipping read mapping step!')

        if args.bwa:
            bamfile = f'{sample_dir}.bam'
            logger.debug(f'bamfile is: {bamfile}')
            if not utils.file_exists_and_not_empty(bamfile):
                fill = textwrap.fill(f'{"[ERROR]:":10} The parameter "--start_from {args.start_from}" has been provided '
                                     f'together with "--bwa", but no existing BWA bamfile with filename "{bamfile}" '
                                     f'can be found. Are you sure you have run the pipeline mapping step already for '
                                     f'this sample?',
                                     width=90, subsequent_indent=' ' * 11)
                logger.info(fill)
                sys.exit(1)
        elif args.blast:
            blastx_outputfile = f'{sample_dir}.blastx'
            logger.debug(f'blastx outputfile is: {blastx_outputfile}')
            if not utils.file_exists_and_not_empty(blastx_outputfile):
                fill = textwrap.fill(f'{"[ERROR]:":10} The parameter "--start_from {args.start_from}" has been provided '
                                     f' and BLASTx/DIAMOND is used for read mapping, but no existing *.blastx output '
                                     f'file with filename "{blastx_outputfile}" can be found. Are you sure you have '
                                     f'run the pipeline mapping step already for this sample?',
                                     width=90, subsequent_indent=' ' * 11)
                logger.info(fill)
                sys.exit(1)
    else:
        if args.bwa:  # map reads to nucleotide targets with BWA
            if args.unpaired:
                # Note that unpaired_readfile is a single path to the file:
                bwa(unpaired_readfile, targetfile, sample_dir, cpu=cpu, unpaired=True, logger=logger)
            # Note that readfiles is a list of one (single-end) or two (paired-end) paths to read files:
            bamfile = bwa(readfiles, targetfile, sample_dir, cpu=cpu, logger=logger)
            if not bamfile:
                logger.error(f'{"[ERROR]:":10} Something went wrong with the BWA step, exiting. Check the log file for '
                             f'sample {sample_dir}!')
                return
            logger.debug(f'bamfile is: {bamfile}')

        elif args.blast:  # map reads to protein targets with BLASTx
            if args.unpaired:
                blastx(unpaired_readfile, targetfile, args.evalue, sample_dir, cpu=cpu,
                       max_target_seqs=args.max_target_seqs, unpaired=True, logger=logger, diamond=args.diamond,
                       diamond_sensitivity=args.diamond_sensitivity)

            blastx_outputfile = blastx(readfiles, targetfile, args.evalue, sample_dir, cpu=cpu,
                                       max_target_seqs=args.max_target_seqs, logger=logger, diamond=args.diamond,
                                       diamond_sensitivity=args.diamond_sensitivity)

            if not blastx_outputfile:
                logger.error(f'{"[ERROR]:":10} Something went wrong with the Blastx step, exiting. Check the log file '
                             f'for sample {sample_dir}!')
                return
        else:
            sys.exit(f'Can not determine whether BWA or BLASTx option is supplied, exiting...')

    if args.end_with in ['map_reads']:
        logger.info(f'{"[INFO]:":10} Parameter "--end_with {args.end_with}" supplied. Exiting assembly pipeline.')
        sys.exit()

    ####################################################################################################################
    # Distribute reads to gene directories from either BLASTx or BWA mapping
    ####################################################################################################################

    if args.start_from not in ['map_reads', 'distribute_reads']:
        logger.info(f'{"[INFO]:":10} Parameter "--start_from {args.start_from}" supplied, skipping read and target '
                    f'distribution step!')
        pre_existing_fastas = glob.glob('./*/*_interleaved.fasta') + glob.glob('./*/*_unpaired.fasta')
        if len(pre_existing_fastas) == 0:
            fill = textwrap.fill(f'{"[ERROR]:":10} The parameter "--start_from {args.start_from}" has been '
                                 f'provided but no distributed reads (*_interleaved.fasta and/or *_unpaired.fasta) can '
                                 f'be found for any gene. Are you sure you have run the pipeline read distribution '
                                 f'step already for this sample?',
                                 width=90, subsequent_indent=' ' * 11)
            logger.info(fill)
            sys.exit(1)
    else:
        pre_existing_fastas = glob.glob('./*/*_interleaved.fasta') + glob.glob('./*/*_unpaired.fasta')
        for fasta in pre_existing_fastas:
            os.remove(fasta)

        if args.bwa:
            distribute_bwa(bamfile, readfiles, targetfile, target, unpaired_readfile, args.exclude,
                           merged=args.merged, low_mem=args.distribute_low_mem, logger=logger)
        else:  # distribute BLASTx results
            distribute_blastx(blastx_outputfile, readfiles, targetfile, target, unpaired_readfile, args.exclude,
                              merged=args.merged, low_mem=args.distribute_low_mem, logger=logger)

    # Note that HybPiper expects either paired-end readfiles (parameter --readfiles) and an optional file of unpaired
    #  reads (parameter --unpaired), or a single file of unpaired reads (parameter --readfiles). For each scenario,
    #  the unpaired readfile is written to an *_unpaired.fasta file.
    if len(readfiles) == 2:
        genes = [x for x in os.listdir('.') if os.path.isfile(os.path.join(x, x + '_interleaved.fasta'))]
    else:
        genes = [x for x in os.listdir('.') if os.path.isfile(os.path.join(x, x + '_unpaired.fasta'))]
    if len(genes) == 0:
        logger.error(f'{"[ERROR]:":10} No genes with reads, exiting!')
        return

    if args.end_with in ['distribute_reads']:
        logger.info(f'{"[INFO]:":10} Parameter "--end_with {args.end_with}" supplied. Exiting assembly pipeline.')
        sys.exit()

    ####################################################################################################################
    # Assemble reads using SPAdes
    ####################################################################################################################

    if args.start_from not in ['map_reads', 'distribute_reads', 'assemble_reads']:
        logger.info(f'{"[INFO]:":10} Parameter "--start_from {args.start_from}" supplied, skipping read assembly step!')
        pre_existing_assemblies = glob.glob('./*/*_contigs.fasta')
        if len(pre_existing_assemblies) == 0:
            fill = textwrap.fill(f'{"[ERROR]:":10} The parameter "--start_from {args.start_from}" has been provided but '
                                 f'no existing SPAdes assembly files (*_contigs.fasta) can be found for any gene. Are '
                                 f'you sure you have run the pipeline read assembly step already for this sample?',
                                 width=90, subsequent_indent=' ' * 11)
            logger.info(fill)
            sys.exit(1)
    else:
        # If the --merged flag is provided, merge reads for SPAdes assembly
        if args.merged:
            logger.info(f'{"[INFO]:":10} Merging reads for SPAdes assembly')
            for gene in genes:
                interleaved_reads_for_merged = f'{gene}/{gene}_interleaved.fastq'
                logger.debug(f'interleaved_reads_for_merged file is {interleaved_reads_for_merged}\n')
                merged_out = f'{gene}/{gene}_merged.fastq'
                unmerged_out = f'{gene}/{gene}_unmerged.fastq'
                bbmerge_command = f'bbmerge.sh interleaved=true in={interleaved_reads_for_merged} out={merged_out}  ' \
                                  f'outu={unmerged_out}'
                try:
                    result = subprocess.run(bbmerge_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                            universal_newlines=True, check=True)
                    logger.debug(f'bbmerge check_returncode() is: {result.check_returncode()}')
                    logger.debug(f'bbmerge paired stdout is: {result.stdout}')
                    logger.debug(f'bbmerge paired stderr is: {result.stderr}')

                except subprocess.CalledProcessError as exc:
                    logger.error(f'bbmerge paired FAILED. Output is: {exc}')
                    logger.error(f'bbmerge paired stdout is: {exc.stdout}')
                    logger.error(f'bbmerge paired stderr is: {exc.stderr}')
                    sys.exit('There was an issue when merging reads. Check read files!')

        if len(readfiles) == 1:
            spades_genelist = spades(genes,
                                     cov_cutoff=args.cov_cutoff,
                                     cpu=cpu,
                                     kvals=args.kvals,
                                     paired=False,
                                     timeout=args.timeout_assemble,
                                     logger=logger,
                                     keep_folder=args.keep_intermediate_files,
                                     single_cell_mode=args.spades_single_cell)
        elif len(readfiles) == 2:
            unpaired = True if unpaired_readfile else False
            spades_genelist = spades(genes,
                                     cov_cutoff=args.cov_cutoff,
                                     cpu=cpu,
                                     kvals=args.kvals,
                                     paired=True,
                                     timeout=args.timeout_assemble,
                                     merged=args.merged,
                                     unpaired=unpaired,
                                     logger=logger,
                                     keep_folder=args.keep_intermediate_files,
                                     single_cell_mode=args.spades_single_cell)

        else:
            logger.error('ERROR: Please specify either one (unpaired) or two (paired) read files! Exiting!')
            return
        if not spades_genelist:
            logger.error('ERROR: No genes had assembled contigs! Exiting!')
            return

    if args.end_with in ['assemble_reads']:
        logger.info(f'{"[INFO]:":10} Parameter "--end_with {args.end_with}" supplied. Exiting assembly pipeline.')
        sys.exit()

    ####################################################################################################################
    # Run Exonerate on the assembled SPAdes contigs, and intronerate() if flag --run_intronerate is supplied:
    ####################################################################################################################

    genes = [x.rstrip() for x in open('exonerate_genelist.txt').readlines()]
    # Remove any pre-existing directories:
    for g in genes:
        if os.path.isdir(os.path.join(g, sample_dir)):
            shutil.rmtree(os.path.join(g, sample_dir))
    if os.path.isfile('genes_with_seqs.txt'):
        os.remove('genes_with_seqs.txt')
    if len(genes) == 0:
        logger.error(f'{"[ERROR]:":10} No genes recovered for {sample_dir}!')
        return 1

    exonerate_multiprocessing(genes,
                              sample_dir,
                              thresh=args.thresh,
                              paralog_warning_min_length_percentage=args.paralog_min_length_percentage,
                              depth_multiplier=args.depth_multiplier,
                              no_stitched_contig=args.no_stitched_contig,
                              stitched_contig_pad_n=args.stitched_contig_pad_n,
                              chimera_check=args.chimera_check,
                              bbmap_memory=args.bbmap_memory,
                              bbmap_subfilter=args.bbmap_subfilter,
                              chimeric_stitched_contig_edit_distance=args.chimeric_stitched_contig_edit_distance,
                              chimeric_stitched_contig_discordant_reads_cutoff=
                              args.chimeric_stitched_contig_discordant_reads_cutoff,
                              bbmap_threads=args.bbmap_threads,
                              pool_threads=cpu,
                              logger=logger,
                              no_intronerate=args.no_intronerate,
                              no_padding_supercontigs=args.no_padding_supercontigs,
                              keep_intermediate_files=args.keep_intermediate_files,
                              exonerate_contigs_timeout=args.timeout_exonerate_contigs,
                              exonerate_hit_sliding_window_size=args.exonerate_hit_sliding_window_size,
                              exonerate_hit_sliding_window_thresh=args.exonerate_hit_sliding_window_thresh,
                              exonerate_skip_frameshifts=args.skip_frameshifts,
                              exonerate_skip_internal_stops=args.skip_internal_stops,
                              exonerate_skip_terminal_stops=args.skip_terminal_stops,
                              verbose_logging=args.verbose_logging)

    ####################################################################################################################
    # Collate all stitched contig and putative chimera read reports
    ####################################################################################################################
    logger.info(f'\n{"[INFO]:":10} Generated sequences from {len(open("genes_with_seqs.txt").readlines())} genes!')
    num_genes_with_stop_codons = len(open(f'{sample_dir}_genes_with_non_terminal_stop_codons.txt').readlines())

    if num_genes_with_stop_codons:
        fill_1 = textwrap.fill(f'{"[WARNING]:":10} {num_genes_with_stop_codons} genes contain internal stop codons. '
                               f'See file "{sample_dir}_genes_with_non_terminal_stop_codons.txt" for a list of gene '
                               f'names, and visit the wiki at the following link to view troubleshooting '
                               f'recommendations:',
                               width=90, subsequent_indent=" " * 11)

        fill_2 = (f'{" " * 10} https://github.com/mossmatters/HybPiper/wiki/Troubleshooting,-common-issues,'
                  f'-and-recommendations#31-sequences-containing-stop-codons')

        logger.warning(f'{fill_1}\n{fill_2}')

    # Stitched contigs:
    collate_stitched_contig_reports = [x for x in glob.glob(f'*/{sample_dir}/genes_with_stitched_contig.csv')]
    with open(f'{sample_dir}_genes_with_stitched_contig.csv', 'w') as genes_with_stitched_contig_handle:
        for report_file in collate_stitched_contig_reports:
            with open(report_file, 'r') as report_handle:
                lines = report_handle.readlines()
                genes_with_stitched_contig_handle.write('\n'.join(lines))

    # Putative chimeras:
    collate_putative_chimeras_reports = [x for x in glob.glob(f'*/{sample_dir}/putative_chimeric_stitched_contig.csv')]
    with open(f'{sample_dir}_genes_derived_from_putative_chimeric_stitched_contig.csv',
              'w') as genes_with_chimeras_handle:
        for report_file in collate_putative_chimeras_reports:
            with open(report_file, 'r') as report_handle:
                lines = report_handle.readlines()
                genes_with_chimeras_handle.write('\n'.join(lines))


    ####################################################################################################################
    # Report paralog warnings and write paralog warning files
    ####################################################################################################################

    # Collate report for long paralogs, and write warning to screen:
    paralog_warnings_long = [x for x in glob.glob(f'*/{sample_dir}/paralog_warning_long.txt')]
    with open(f'{sample_dir}_genes_with_long_paralog_warnings.txt', 'w') as long_paralogs_handle:
        for warning_file in paralog_warnings_long:
            with open(warning_file, 'r') as warning_handle:
                report_line = warning_handle.readline().split()[0]  # only recover gene name
                long_paralogs_handle.write(f'{report_line}\n')
    logger.info(f'{"[WARNING]:":10} Potential long paralogs detected for {len(paralog_warnings_long)} genes!')

    # Collate report for paralogs via SPAdes contig depth, and write warning to screen:
    paralog_warnings_short = [x for x in glob.glob(f'*/{sample_dir}/paralog_warning_by_contig_depth.txt')]
    paralog_warnings_short_true = 0
    with open(f'{sample_dir}_genes_with_paralog_warnings_by_contig_depth.csv', 'w') as depth_paralogs_handle:
        for warning_file in paralog_warnings_short:
            with open(warning_file, 'r') as warning_handle:
                report_line = warning_handle.readline()
                if report_line.split()[-1] == 'True':
                    paralog_warnings_short_true += 1
                depth_paralogs_handle.write(report_line)
    logger.info(f'{"[WARNING]:":10} Potential paralogs detected via contig depth for'
                f' {paralog_warnings_short_true} genes!')

    logger.info(f'\nFinished running "hybpiper assemble" for sample {sample_dir}!\n')
