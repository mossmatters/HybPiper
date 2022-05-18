#!/usr/bin/env python

"""
HybPiper Version 2.0 release candidate (May 2022)

########################################################################################################################
############################################## NOTES ON VERSION 2.0 ####################################################
########################################################################################################################

After installation of the pipeline, all pipeline commands are now accessed via the main command 'hybpiper',
followed by a subcommand to run different parts of the pipeline. The available subcommands can be viewed by typing
'hybpiper -h' or 'hybpiper --help'. They are:

    assemble            Assemble gene, intron, and supercontig sequences
    stats               Gather statistics about the HybPiper run(s)
    retrieve_sequences  Retrieve sequences generated from multiple runs of HybPiper
    recovery_heatmap    Create a gene recovery heatmap for the HybPiper run
    paralog_retriever   Retrieve paralog sequences for a given gene, for all samples
    check_dependencies  Run a check for all pipeline dependencies and exit
    check_targetfile    Check the target file for sequences with low-complexity regions, then exit

To view available parameters and help for any subcommand, simply type e.g. 'hybpiper assemble -h'.

==> NOTE <==
The script 'read_first.py' no longer exists, and has been replaced by the subcommand 'assemble'. So,
if you had previously run 'reads_first.py' on a sample using the command e.g.:

    python /<path_to>/reads_first.py -b test_targets.fasta -r NZ281_R*_test.fastq --prefix NZ281 --bwa

...this is now replaced by the command:

    hybpiper assemble -t_dna test_targets.fasta -r NZ281_R*_test.fastq --prefix NZ281 --bwa

==> NOTE <==
The recovery of introns and supercontigs, previously achieved via the script 'intronerate.py',
is now incorporated in to the 'hybpiper assemble' command. It can be enabled using the flag
'--run_intronerate', e.g.:

    hybpiper assemble -t_dna test_targets.fasta -r NZ281_R*_test.fastq --prefix NZ281 --bwa --run_intronerate

==> NOTE <==
The command/script 'get_seq_lengths.py' no longer exists, and this functionality has been incorporated in to
the command 'hybpiper stats'. The sequence length details that were previously printed to screen are now written to
the file 'seq_lengths.tsv', by default. Similarly, the stats details that were previously written to screen by
'hybpiper_stats.py' are now written to the file 'hybpiper_stats.tsv', by default.

For full details of all commands and changes, please read the Wiki page at
https://github.com/mossmatters/HybPiper/wiki and the change log at
https://github.com/mossmatters/HybPiper/blob/master/change_log.md.

########################################################################################################################

"""

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
import datetime
import multiprocessing
from multiprocessing import Manager
from concurrent.futures import wait, as_completed, TimeoutError, CancelledError
import pkg_resources
import time
import signal

# f-strings will produce a 'SyntaxError: invalid syntax' error if not supported by Python version:
f'HybPiper requires Python 3.6 or higher.'

# Import non-standard-library modules:
unsuccessful_imports = []
try:
    from Bio import SeqIO, SeqRecord
except ImportError:
    unsuccessful_imports.append('Bio')
try:
    import progressbar
except ImportError:
    unsuccessful_imports.append('progressbar')
try:
    import seaborn
except ImportError:
    unsuccessful_imports.append('seaborn')
try:
    import matplotlib
except ImportError:
    unsuccessful_imports.append('matplotlib')
try:
    import pandas
except ImportError:
    unsuccessful_imports.append('pandas')
try:
    import pebble
except ImportError:
    unsuccessful_imports.append('pebble')
try:
    import scipy
except ImportError:
    unsuccessful_imports.append('scipy')
try:
    import psutil
except ImportError:
    unsuccessful_imports.append('psutil')

if unsuccessful_imports:
    package_list = '\n'.join(unsuccessful_imports)
    sys.exit(f'The required Python packages are not found:\n\n{package_list}\n\nAre they installed for the Python '
             f'installation used to run HybPiper?')

# Check that user has the minimum required version of Biopython (1.80):
biopython_version_print = pkg_resources.get_distribution('biopython').version
biopython_version = [int(value) for value in re.split('[.]', biopython_version_print)[:2]]
if biopython_version[0:2] < [1, 80]:
    sys.exit(f'HybPiper required Biopython version 1.80 or above. You are using version {biopython_version_print}. '
             f'Please update your Biopython for the Python installation used to run HybPiper!')

# Import HybPiper modules:
from hybpiper import distribute_reads_to_targets
from hybpiper import distribute_targets
from hybpiper import spades_runner
from hybpiper import exonerate_hits
from hybpiper import hybpiper_stats
from hybpiper import retrieve_sequences
from hybpiper import paralog_retriever
from hybpiper import gene_recovery_heatmap
from hybpiper import hybpiper_subparsers
from hybpiper import utils
from hybpiper.utils import log_or_print


########################################################################################################################
# Define functions
########################################################################################################################

def setup_logger(name, log_file, console_level=logging.INFO, file_level=logging.DEBUG,
                 logger_object_level=logging.DEBUG):
    """
    Function to create a logger instance.

    By default, logs level DEBUG and above to file.
    By default, logs level INFO and above to stderr and file.

    :param str name: name for the logger instance
    :param str log_file: filename for log file
    :param str console_level: logger level for logging to console
    :param str file_level: logger level for logging to file
    :param str logger_object_level: logger level for logger object
    :return: logging.Logger: logger object
    """

    # Get date and time string for log filename:
    date_and_time = datetime.datetime.now().strftime("%Y-%m-%d-%H_%M_%S")

    # Log to file:
    file_handler = logging.FileHandler(f'{log_file}_{date_and_time}.log', mode='w')
    file_handler.setLevel(file_level)
    file_format = logging.Formatter('%(asctime)s - %(filename)s - %(name)s - %(funcName)s - %(levelname)s - %('
                                    'message)s')
    file_handler.setFormatter(file_format)

    # Log to Terminal (stderr):
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(console_level)
    console_format = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_format)

    # Setup logger:
    logger_object = logging.getLogger(name)
    logger_object.setLevel(logger_object_level)  # Default level is 'WARNING'

    # Add handlers to the logger
    logger_object.addHandler(console_handler)
    logger_object.addHandler(file_handler)

    return logger_object


def check_target_file_headers_and_duplicate_names(targetfile, logger=None):
    """
    - Checks target-file fasta header formatting ("taxon*-unique_gene_ID").
    - Checks for duplicate gene names in the targetfile.
    - Reports the number of unique genes (each can have multiple representatives) in the targetfile.

    :param str targetfile: path to the targetfile
    :param logging.Logger logger: a logger object
    :return:
    """

    log_or_print(f'{"[INFO]:":10} Checking target file FASTA header formatting...', logger=logger)

    gene_lists = defaultdict(list)
    with open(targetfile, 'r') as target_file_handle:
        seqs = list(SeqIO.parse(target_file_handle, 'fasta'))
        incorrectly_formatted_fasta_headers = []
        check_for_duplicate_genes_dict = {}
        for seq in seqs:
            if seq.name in check_for_duplicate_genes_dict:
                check_for_duplicate_genes_dict[seq.name] += 1
            else:
                check_for_duplicate_genes_dict[seq.name] = 1
            if not re.match('.+-[^-]+', seq.name):
                incorrectly_formatted_fasta_headers.append(seq.name)
            gene_id = re.split('-', seq.name)[-1]
            gene_lists[gene_id].append(seq)

    if incorrectly_formatted_fasta_headers:
        seq_list = ' '.join(incorrectly_formatted_fasta_headers)
        log_or_print(f'{"[ERROR!]:":10} The following sequences in your target file have incorrectly formatted fasta '
                     f'headers:\n', logger=logger, logger_level='error')
        fill = textwrap.fill(f'{seq_list}')
        log_or_print(textwrap.indent(fill, ' ' * 11), logger=logger)
        log_or_print('', logger=logger)
        sys.exit(1)  # target file fasta header formatting should be fixed!
    else:
        log_or_print(f'{"[INFO]:":10} The target file FASTA header formatting looks good!', logger=logger)

    # Check for duplicated gene names:
    duplicated_genes = []
    for gene, gene_count in check_for_duplicate_genes_dict.items():
        if gene_count > 1:
            duplicated_genes.append(gene)
    if duplicated_genes:
        gene_list = ' '.join(duplicated_genes)
        log_or_print(f'{"[ERROR!]:":10} The following sequences in your target file occur more than once:\n',
                     logger=logger, logger_level='error')
        fill = textwrap.fill(f'{gene_list}')
        log_or_print(textwrap.indent(fill, ' ' * 11), logger=logger)
        log_or_print(f'\nPlease remove duplicate genes before running HybPiper!', logger=logger, logger_level='error')
        sys.exit(1)  # duplicate genes in target file should be removed!

    # Report the number of unique genes represented in the target file:
    log_or_print(f'{"[INFO]:":10} The target file contains at least one sequence for {len(gene_lists)} '
                 f'unique genes.', logger=logger)


def check_target_file_stop_codons_and_multiple_of_three(targetfile, translate_target_file=False, logger=None):
    """

    :param targetfile:
    :param translate_target_file:
    :param logger:
    :return:
    """

    with open(targetfile, 'r') as target_file_handle:
        seqs = list(SeqIO.parse(target_file_handle, 'fasta'))

    translated_seqs_to_write = []
    seqs_needed_padding_dict = defaultdict(list)
    seqs_with_stop_codons_dict = defaultdict(list)
    seqs_with_stop_codon_last_five_amino_acids_dict = defaultdict(list)

    if translate_target_file:
        for seq in seqs:
            gene_name = seq.name.split('-')[-1]
            sequence, needed_padding = utils.pad_seq(seq)
            translated_seq = sequence.seq.translate()
            if translate_target_file:
                record = SeqRecord.SeqRecord(translated_seq, id=seq.id, description='')
                translated_seqs_to_write.append(record)
            num_stop_codons = translated_seq.count('*')

            if needed_padding:
                seqs_needed_padding_dict[gene_name].append(seq)

            if num_stop_codons == 1 and re.search('[*]', str(translated_seq)[-5:]):
                seqs_with_stop_codon_last_five_amino_acids_dict[gene_name].append(seq)
            elif num_stop_codons >= 1:
                seqs_with_stop_codons_dict[gene_name].append(seq)

        if seqs_with_stop_codon_last_five_amino_acids_dict:
            seq_list = [seq.name for gene_name, target_file_sequence_list in
                        seqs_with_stop_codon_last_five_amino_acids_dict.items() for seq in target_file_sequence_list]
            fill = textwrap.fill(
                f'{"[INFO]:":10} There are {len(seq_list)} sequences in your target file that contain a single stop '
                f'codon in the last five amino-acids of the C-terminal end. These are: ', width=90,
                subsequent_indent=' ' * 11)
            log_or_print(f'{fill}\n', logger=logger, logger_level='debug')
            fill = textwrap.fill(f'{", ".join(seq_list)}', width=90, initial_indent=' ' * 11,
                                 subsequent_indent=' ' * 11, break_on_hyphens=False)
            log_or_print(f'{fill}\n', logger=logger, logger_level='debug')

        if seqs_with_stop_codons_dict:
            seq_list = [seq.name for gene_name, target_file_sequence_list in seqs_with_stop_codons_dict.items() for seq
                        in target_file_sequence_list]
            fill = textwrap.fill(
                f'{"[WARNING]:":10} There are {len(seq_list)} sequences in your target file that contain unexpected '
                f'stop codons when translated in the first forwards frame. If your target file contains only '
                f'protein-coding sequences, please check these sequences. Sequence names can be found in the sample '
                f'log file (if running "hybpiper assemble") or printed below (if running "hybpiper check_targetfile").',
                width=90, subsequent_indent=' ' * 11)
            log_or_print(f'{fill}\n', logger=logger)
            fill = textwrap.fill(f'{", ".join(seq_list)}', width=90, initial_indent=' ' * 11,
                                 subsequent_indent=' ' * 11, break_on_hyphens=False)
            log_or_print(f'{fill}\n', logger=logger)

        if seqs_needed_padding_dict:
            seq_list = [seq.name for gene_name, target_file_sequence_list in seqs_needed_padding_dict.items() for seq
                        in target_file_sequence_list]
            fill = textwrap.fill(
                f'{"[WARNING]:":10} There are {len(seq_list)} sequences in your target file that are not multiples of '
                f'three. If your target file contains only protein-coding sequences, please check these sequences. '
                f'Sequence names can be found in the sample log file (if running "hybpiper assemble") or printed '
                f'below (if running "hybpiper check_targetfile")..', width=90, subsequent_indent=' ' * 11)
            log_or_print(f'{fill}\n', logger=logger)
            fill = textwrap.fill(f'{", ".join(seq_list)}', width=90, initial_indent=' ' * 11,
                                 subsequent_indent=' ' * 11, break_on_hyphens=False)
            log_or_print(f'{fill}\n', logger=logger)

    return translated_seqs_to_write


def check_targetfile(targetfile, targetfile_type, using_bwa, logger=None):
    """
    - Checks target-file fasta header formatting ("taxon*-unique_gene_ID").
    - Checks for duplicate gene names in the targetfile.
    - Reports the number of unique genes (each can have multiple representatives) in the targetfile.
    - Checks that seqs in target file can be translated from the first codon position in the forwards frame (multiple of
      three, no unexpected stop codons), and logs a warning if not.
    - If targetfile is DNA but using_bwa is False, translate the targetfile and return the path

    :param str targetfile: path to the targetfile
    :param str targetfile_type: string describing target file sequence type i.e 'DNA' or 'protein'
    :param bool using_bwa: True if the --bwa flag is used; a nucleotide target file is expected in this case
    :param logging.Logger logger: a logger object
    :return: None, str: NoneType or path to the translated targetfile
    """

    target_file_path, target_file_name = os.path.split(targetfile)
    file_name, ext = os.path.splitext(target_file_name)

    # Check target file header and duplicate gene names:
    check_target_file_headers_and_duplicate_names(targetfile, logger=logger)

    # Detect whether the target file is DNA or amino-acid:
    translate_target_file = False
    if using_bwa and targetfile_type == 'protein':
        sys.exit(f'{"[ERROR]:":10} You have specified that your target file contains protein sequences but provided '
                 f'the flag --bwa. You need a nucleotide target file to use BWA for read mapping!')
    elif not using_bwa and targetfile_type == 'DNA':
        fill = textwrap.fill(f'{"[WARNING]:":10} You have specified that your target file contains DNA sequences, '
                             f'but BLASTx or DIAMOND has been selected for read mapping. Translating the target '
                             f'file...', width=90, subsequent_indent=' ' * 11)
        logger.info(f'{fill}')
        translate_target_file = True

    # Check that seqs in target file can be translated from the first codon position in the forwards frame:
    if using_bwa or translate_target_file:  # i.e. it's not a protein file
        translated_seqs_to_write = check_target_file_stop_codons_and_multiple_of_three(
            targetfile,
            translate_target_file=True,  # Set manually
            logger=logger)

        if translate_target_file:
            translated_target_file = f'{target_file_path}/{file_name}_translated{ext}'
            fill = utils.fill_forward_slash(f'{"[INFO]:":10} Writing a translated target file to:'
                                            f' {translated_target_file}', width=90, subsequent_indent=' ' * 11,
                                            break_long_words=False, break_on_forward_slash=True)
            logger.info(f'{fill}')

            with open(f'{translated_target_file}', 'w') as translated_handle:
                SeqIO.write(translated_seqs_to_write, translated_handle, 'fasta')
            targetfile = translated_target_file  # i.e. use translated file for return value

    return targetfile


def bwa(readfiles, targetfile, basename, cpu, unpaired=False, logger=None):
    """
    Conduct a BWA search of input reads against the targetfile.

    :param str/list readfiles: list one or more read files used as input to the pipeline, or path to unpaired read file
    :param str targetfile: path to targetfile (i.e. the target file)
    :param str basename: directory name for sample
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
            fill = textwrap.fill(f'{"[CMD]:":10} {make_bwa_index_cmd}', width=90, subsequent_indent=' ' * 11,
                                 break_long_words=False, break_on_hyphens=False)
            logger.info(f'{fill}')

            try:
                result = subprocess.run(make_bwa_index_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True, check=True)
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

    bwa_commands = ['time bwa mem', '-t', str(cpu), db_file, bwa_fastq, ' | samtools view -h -b -S - > ']
    if unpaired:
        bwa_commands.append(f'{basename}_unpaired.bam')
    else:
        bwa_commands.append(f'{basename}.bam')
    full_command = ' '.join(bwa_commands)
    fill = utils.fill_forward_slash(f'{"[CMD]:":10} {full_command}', width=90, subsequent_indent=' ' * 11,
                                    break_long_words=False, break_on_forward_slash=True)
    logger.info(f'{fill}')

    try:
        result = subprocess.run(full_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True, check=True)
        logger.debug(f'BWA mapping check_returncode() is: {result.check_returncode()}')
        logger.debug(f'BWA mapping stdout is: {result.stdout}')
        logger.debug(f'BWA mapping stderr is: {result.stderr}')

    except subprocess.CalledProcessError as exc:
        logger.error(f'BWA mapping FAILED. Output is: {exc}')
        logger.error(f'BWA mapping stdout is: {exc.stdout}')
        logger.error(f'BWA mapping stderr is: {exc.stderr}')
        return None

    return f'{basename}.bam'  # No return for {basename}_unpaired.bam?


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

    targetfile_basename = os.path.basename(targetfile)

    if os.path.isfile(targetfile):
        if os.path.isfile(f'{targetfile_basename}.psq'):
            db_file = targetfile_basename
            logger.debug(f'Using existing BLAST database. db_file is: {db_file}')
        elif os.path.isfile(f'{targetfile_basename}.diamond'):
            db_file = targetfile_basename
            logger.debug(f'Using existing DIAMOND BLAST database. db_file is: {db_file}')
        else:
            logger.info(f'{"[INFO]:":10} Making protein blastdb in current directory.')
            if os.path.split(targetfile)[0]:
                shutil.copy(targetfile, '.')
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
                result = subprocess.run(makeblastdb_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True, check=True)
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
        if file_extension == '.gz':
            logger.debug(f'Processing gzipped file {os.path.basename(read_file)}')
            pipe_cmd = f"gunzip -c {read_file} | awk '{{if(NR % 4 == 1 || NR % 4 == 2) {{sub(/@/, \">\"); print; }} " \
                       f"}}'"
        else:
            pipe_cmd = f"cat {read_file} | awk '{{if(NR % 4 == 1 || NR % 4 == 2) {{sub(/@/, \">\"); print; }} }}'"
        if diamond and diamond_sensitivity:
            blastx_command = f'diamond blastx --db {db_file} --query - --evalue {evalue} --outfmt 6 --max-target-seqs' \
                             f' {max_target_seqs} --{diamond_sensitivity}'
        elif diamond:
            blastx_command = f'diamond blastx --db {db_file} --query - --evalue {evalue} --outfmt 6 --max-target-seqs' \
                             f' {max_target_seqs}'
        else:
            blastx_command = f'blastx -db {db_file} -query - -evalue {evalue} -outfmt 6 -max_target_seqs' \
                             f' {max_target_seqs}'

        full_command = f"time {pipe_cmd} | parallel -j {cpu} -k --block 200K --recstart '>' --pipe '{blastx_command}' >>" \
                       f" {basename}_unpaired.blastx"

        fill = utils.fill_forward_slash(f'{"[CMD]:":10} {full_command}', width=90, subsequent_indent=' ' * 11,
                                        break_long_words=False, break_on_forward_slash=True)
        logger.info(f'{fill}')

        try:
            result = subprocess.run(full_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True, check=True)
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
            if file_extension == '.gz':
                logger.debug(f'Processing gzipped file {os.path.basename(read_file)}')
                pipe_cmd = f"gunzip -c {read_file} | awk '{{if(NR % 4 == 1 || NR % 4 == 2) {{sub(/@/, \">\"); print; " \
                           f"}} }}'"
            else:
                pipe_cmd = f"cat {read_file} | awk '{{if(NR % 4 == 1 || NR % 4 == 2) {{sub(/@/, \">\"); print; }} }}'"

            if diamond and diamond_sensitivity:
                blastx_command = f'diamond blastx --db {db_file} --query - --evalue {evalue} --outfmt 6 ' \
                                 f'--max-target-seqs {max_target_seqs} --{diamond_sensitivity}'
            elif diamond:
                blastx_command = f'diamond blastx --db {db_file} --query - --evalue {evalue} --outfmt 6 ' \
                                 f'--max-target-seqs {max_target_seqs}'
            else:
                blastx_command = f'blastx -db {db_file} -query - -evalue {evalue} -outfmt 6 -max_target_seqs' \
                                 f' {max_target_seqs}'

            full_command = f"time {pipe_cmd} | parallel -j {cpu} -k --block 200K --recstart '>' --pipe " \
                           f"'{blastx_command}' >> {basename}.blastx"

            fill = utils.fill_forward_slash(f'{"[CMD]:":10} {full_command}', width=90, subsequent_indent=' ' * 11,
                                            break_long_words=False, break_on_forward_slash=True)
            logger.info(f'{fill}')

            try:
                result = subprocess.run(full_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True, check=True)
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
                      merged=False, hi_mem=False, logger=None):
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
    :param bool hi_mem: If True, reads to distribute will be saved in a dictionary and written once; used more RAM
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

    distribute_reads_to_targets.distribute_reads(readfiles, read_hit_dict, merged=merged, single_end=single_end,
                                                 hi_mem=hi_mem)

    if unpaired_readfile:
        up_blastx_outputfile = blastx_outputfile.replace('.blastx', '_unpaired.blastx')
        read_hit_dict_unpaired = distribute_reads_to_targets.read_sorting_blastx(up_blastx_outputfile)
        logger.info(f'{"[INFO]:":10} In total, {len(read_hit_dict_unpaired)} reads from the unpaired read file will be '
                    f'distributed to gene directories')
        distribute_reads_to_targets.distribute_reads([unpaired_readfile], read_hit_dict_unpaired,
                                                     unpaired_readfile=unpaired_readfile, hi_mem=hi_mem)

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

    besthits = distribute_targets.tailored_target_blast(blastx_outputfile, unpaired_bool, exclude_string)
    distribute_targets.distribute_targets(targetfile, delim='-', besthits=besthits, translate=False,
                                          target=target_string)
    return None


def distribute_bwa(bamfile, readfiles, targetfile, target=None, unpaired_readfile=None, exclude=None, merged=False,
                   hi_mem=False, logger=None):
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
    :param bool hi_mem: If True, reads to distribute will be saved in a dictionary and written once; used more RAM
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

    distribute_reads_to_targets.distribute_reads(readfiles, read_hit_dict, merged=merged, single_end=single_end,
                                                 hi_mem=hi_mem)

    if unpaired_readfile:
        up_bamfile = bamfile.replace('.bam', '_unpaired.bam')
        read_hit_dict_unpaired = distribute_reads_to_targets.read_sorting_bwa(up_bamfile)
        logger.info(f'{"[INFO]:":10} In total, {len(read_hit_dict_unpaired)} reads from the unpaired read file will be '
                    f'distributed to gene directories')
        distribute_reads_to_targets.distribute_reads([unpaired_readfile], read_hit_dict_unpaired,
                                                     unpaired_readfile=unpaired_readfile, hi_mem=hi_mem)

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

    besthits = distribute_targets.tailored_target_bwa(bamfile, unpaired_bool, exclude_string)
    distribute_targets.distribute_targets(targetfile, delim='-', besthits=besthits, translate=True,
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

    spades_failed = spades_runner.spades_initial('spades_genelist.txt', cov_cutoff=cov_cutoff, cpu=cpu,
                                                 kvals=kvals, paired=paired, timeout=timeout, unpaired=unpaired,
                                                 merged=merged, single_cell_mode=single_cell_mode)
    logger.info(f'{"[INFO]:":10} Finished running initial SPAdes assemblies for all genes with reads!')
    if len(spades_failed) > 0:
        with open('failed_spades.txt', 'w') as failed_spadefile:
            failed_spadefile.write('\n'.join(spades_failed))

        spades_duds = spades_runner.rerun_spades('failed_spades.txt', cov_cutoff=cov_cutoff, cpu=cpu)
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
              basename,
              pid_list,
              thresh=55,
              paralog_warning_min_length_percentage=0.75,
              depth_multiplier=10,
              no_stitched_contig=False,
              bbmap_memory=250,
              bbmap_subfilter=7,
              bbmap_threads=2,
              chimeric_stitched_contig_edit_distance=5,
              chimeric_stitched_contig_discordant_reads_cutoff=5,
              worker_configurer_func=None,
              counter=None,
              lock=None,
              genes_to_process=0,
              intronerate=False,
              no_padding_supercontigs=False,
              keep_intermediate_files=False,
              verbose_logging=False):
    """
    :param str gene_name: name of a gene that had at least one SPAdes contig
    :param str basename: directory name for sample
    :param multiprocessing.managers.ListProxy pid_list: list shared by processes for capturing parent PIDs
    :param int thresh: percent identity threshold for stitching together Exonerate results
    :param float paralog_warning_min_length_percentage: min % of a contig vs ref protein length for a paralog warning
    :param int depth_multiplier: assign long paralog as main if coverage depth <depth_multiplier> other paralogs
    :param bool no_stitched_contig: if True, don't create stitched contigs and just use longest Exonerate hit
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
    :param bool intronerate: if True, run intronerate
    :param bool no_padding_supercontigs: if True, don't pad contig joins in supercontigs with stretches if 10 Ns
    :param bool keep_intermediate_files: if True, keep intermediate files from stitched contig and intronerate()
    processing
    :param bool verbose_logging: if True, log additional information to file
    :return: str gene_name, str prot_length OR None, None
    """

    # get parent PID and add to shared list; this is used to kill child processes on user interrupt:
    pid_list.append(os.getpid())

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
    prefix = exonerate_hits.create_output_directories(f'{gene_name}/{basename}', f'{gene_name}/'
                                                                                 f'{gene_name}_contigs.fasta')
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
            sys.stderr.write(f'\r{"[INFO]:":10} Finished running Exonerate for gene {gene_name}, {counter.value}'
                             f'/{genes_to_process}')

        end = time.time()
        proc_run_time = end - start
        return gene_name, None, proc_run_time  # return gene_name to that log can be re-logged to main log file

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
            discordant_cutoff=chimeric_stitched_contig_discordant_reads_cutoff,
            edit_distance=chimeric_stitched_contig_edit_distance,
            bbmap_subfilter=bbmap_subfilter,
            bbmap_memory=bbmap_memory,
            bbmap_threads=bbmap_threads,
            interleaved_fasta_file=path_to_interleaved_fasta,
            no_stitched_contig=no_stitched_contig,
            spades_assembly_dict=spades_assembly_dict,
            depth_multiplier=depth_multiplier,
            keep_intermediate_files=keep_intermediate_files,
            verbose_logging=verbose_logging)

        if intronerate and exonerate_result and exonerate_result.hits_filtered_by_pct_similarity_dict:
            if verbose_logging:
                logger.debug(f'exonerate_result.hits_subsumed_hits_removed_overlaps_trimmed_dict for gene {gene_name} '
                             f'is: {exonerate_result.hits_subsumed_hits_removed_overlaps_trimmed_dict}')
            logger.debug(f'Running intronerate')
            exonerate_hits.intronerate(exonerate_result,
                                       spades_assembly_dict,
                                       logger=logger,
                                       no_padding_supercontigs=no_padding_supercontigs,
                                       keep_intermediate_files=keep_intermediate_files,
                                       verbose_logging=verbose_logging)
    else:
        exonerate_result = False

    with lock:
        counter.value += 1
        sys.stderr.write(f'\r{"[INFO]:":10} Finished running Exonerate for gene {gene_name}, {counter.value}'
                         f'/{genes_to_process}')

    if not exonerate_text_output or not exonerate_result or not exonerate_result.stitched_contig_seqrecord:
        end = time.time()
        proc_run_time = end - start
        # return gene_name to that exonerate_hits.py log can be re-logged to main log file:
        return gene_name, None, proc_run_time

    end = time.time()
    proc_run_time = end - start

    return gene_name, len(exonerate_result.stitched_contig_seqrecord), proc_run_time


def exonerate_multiprocessing(genes,
                              basename,
                              thresh=55,
                              paralog_warning_min_length_percentage=0.75,
                              pool_threads=None,
                              depth_multiplier=10,
                              no_stitched_contig=False,
                              bbmap_memory=250,
                              bbmap_subfilter=7,
                              bbmap_threads=2,
                              chimeric_stitched_contig_edit_distance=5,
                              chimeric_stitched_contig_discordant_reads_cutoff=5,
                              logger=None,
                              intronerate=False,
                              no_padding_supercontigs=False,
                              keep_intermediate_files=False,
                              exonerate_contigs_timeout=None,
                              verbose_logging=False):
    """
    Runs the function exonerate() using multiprocessing.

    :param list genes: list of genes that had successful SPAdes runs
    :param str basename: directory name for sample
    :param int thresh: percent identity threshold for stitching together Exonerate results
    :param float paralog_warning_min_length_percentage: min % of a contig vs ref protein length for a paralog warning
    :param int pool_threads: number of threads/cpus to use for the pebble.ProcessPool pool
    :param int depth_multiplier: assign long paralog as main if coverage depth <depth_multiplier> other paralogs
    :param bool no_stitched_contig: if True, don't create stitched contig and just use longest Exonerate hit
    :param int bbmap_memory: MB memory (RAM ) to use for bbmap.sh
    :param int bbmap_subfilter: ban alignments with more than this many substitutions
    :param int bbmap_threads: number of threads to use for BBmap when searching for chimeric stitched contigs
    :param int chimeric_stitched_contig_edit_distance: min num differences for a read pair to be flagged as discordant
    :param int chimeric_stitched_contig_discordant_reads_cutoff: min num discordant reads pairs to flag a stitched
    contig as chimeric
    :param logging.Logger logger: a logger object
    :param bool intronerate: if True, intronerate will be run (if a gene is constructed from hits with introns)
    :param bool no_padding_supercontigs: if True, don't pad contig joins in supercontigs with stretches if 10 Ns
    :param bool keep_intermediate_files: if True, keep individual Exonerate logs rather than deleting them after
    re-logging to the main sample log file
    :param int exonerate_contigs_timeout: number of second for pebble.ProcessPool pool.schedule timeout
    :param bool verbose_logging: if True, log additional information to file
    :return:
    """

    logger.debug(f'exonerate_contigs_timeout is: {exonerate_contigs_timeout}')

    logger.info(f'{"[INFO]:":10} Running exonerate_hits for {len(genes)} genes...')
    genes_to_process = len(genes)

    logger.debug(f'exonerate_multiprocessing pool_threads is: {pool_threads}')

    try:
        with pebble.ProcessPool(max_workers=pool_threads) as pool:
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
                                   "intronerate": intronerate,
                                   "no_padding_supercontigs": no_padding_supercontigs,
                                   "keep_intermediate_files": keep_intermediate_files,
                                   "verbose_logging": verbose_logging}

            for gene_name in genes:  # schedule jobs and store each future in a future : gene_name dict
                exonerate_job = pool.schedule(exonerate, args=[gene_name, basename, pid_list],
                                              kwargs=kwargs_for_schedule, timeout=exonerate_contigs_timeout)
                future_results_dict[exonerate_job] = gene_name

            futures_list = [future for future in future_results_dict.keys()]

            # As per-gene Exonerate runs complete, read the gene log, log it to the main logger, delete gene log:
            with open('genes_with_seqs.txt', 'w') as genes_with_seqs_handle:
                for future in as_completed(futures_list):
                    try:
                        gene_name, prot_length, run_time = future.result()
                        if gene_name:  # i.e. log the Exonerate run regardles   s of success
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

                    except TimeoutError as err:
                        logger.debug(f'\nProcess timeout - exonerate() for gene {future_results_dict[future]} took '
                                     f'more than {err.args[1]} seconds to complete and was cancelled')
                        genes_cancelled_due_to_timeout.append(future_results_dict[future])
                    except CancelledError:
                        logger.debug(f'CancelledError raised for gene {future_results_dict[future]}')
                    except Exception as error:
                        genes_cancelled_due_to_errors.append(future_results_dict[future])
                        print(f'For gene {future_results_dict[future]} exonerate() raised: {error}')
                        print(f'error.traceback is: {error.traceback}')  # traceback of the function

        wait(futures_list, return_when="ALL_COMPLETED")  # redundant, but...

        if genes_cancelled_due_to_errors:
            fill = textwrap.fill(f'{"[INFO]:":10} The exonerate_contigs step of the pipeline failed for the '
                                 f'following genes:\n', width=90, subsequent_indent=" " * 11)
            logger.info(fill)
            for gene in genes_cancelled_due_to_errors:
                logger.info(f'{" " * 11}{gene}')

        if genes_cancelled_due_to_timeout:
            fill = textwrap.fill(f'{"[INFO]:":10} The exonerate_contigs step of the pipeline was cancelled for the '
                                 f'following genes, due to exceeding the timeout limit of {exonerate_contigs_timeout} '
                                 f'seconds\n:', width=90, subsequent_indent=" " * 11)
            logger.info(fill)
            for gene in genes_cancelled_due_to_timeout:
                logger.info(f'{" " * 11}{gene}')

            fill = textwrap.fill(f'{"[INFO]:":10} This is most likely caused by many low-complexity reads mapping to '
                                 f'the corresponding gene sequences in the target file, resulting in SPAdes assembly '
                                 f'many (i.e. hundreds) repetitive and low-complexity contigs. Subsequently, '
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


def assemble(args):
    """
    Assemble gene, intron, supercontig and paralog sequences via assemble.py

    :param argparse.Namespace args: argparse namespace with subparser options for function assemble()
    :return None: no return value specified; default is None
    """

    # Get a list of read files from args.readfiles (doesn't include any read file passed in via --unpaired flag):
    readfiles = [os.path.abspath(x) for x in args.readfiles]

    if len(readfiles) > 2:
        sys.exit(f'{"[ERROR]:":10} Please provide a maximum of two read files (R1 and R2) to the -r / '
                 f'--readfiles parameter')

    # Generate a directory for the sample:
    basedir, basename = utils.make_basename(args.readfiles, prefix=args.prefix)

    # Create logger:
    if args.prefix:
        logger = setup_logger(__name__, f'{basename}/{args.prefix}_hybpiper_assemble')
    else:
        logger = setup_logger(__name__, f'{basename}/{os.path.split(readfiles[0])[1].split("_")[0]}_hybpiper_assemble')

    logger.info(f'{"[INFO]:":10} HybPiper was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.info(f'{fill}\n')

    # Get number of cpus/threads for pipeline:
    if args.cpu:
        cpu = args.cpu
        logger.info(f'{"[INFO]:":10} Using {cpu} cpus/threads.')
    else:
        cpu = multiprocessing.cpu_count()  # i.e. use all cpus.
        logger.info(f'{"[INFO]:":10} Number of cpus/threads not specified, using all available ({cpu}).')

    logger.debug(f'args.start_from is: {args.start_from}')

    ####################################################################################################################
    # Check dependencies
    ####################################################################################################################
    if utils.check_dependencies(logger=logger):
        logger.info(f'{"[INFO]:":10} Everything looks good!')
    else:
        logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')

    ####################################################################################################################
    # Check read and target files
    ####################################################################################################################
    # Set target file type and path, and check it exists and isn't empty::
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
    # arg.bwa is false, translate and return the path to translated file:
    targetfile = check_targetfile(targetfile,
                                  targetfile_type,
                                  args.bwa,
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

    # Move in to the sample directory:
    os.chdir(os.path.join(basedir, basename))

    ####################################################################################################################
    # Map reads to target file sequences
    ####################################################################################################################

    if args.start_from not in ['map_reads']:  # 'map_reads' is the default
        logger.info(f'{"[INFO]:":10} Parameter "--start_from {args.start_from}" supplied, skipping read mapping step!')

        if args.bwa:
            bamfile = f'{basename}.bam'
            logger.debug(f'bamfile is: {bamfile}')
            if not utils.file_exists_and_not_empty(bamfile):
                fill = textwrap.fill(f'{"[ERROR]:":10} The parameter "--start_from {args.start_from}" has been provided '
                                     f'together with "--bwa", but no existing BWA bamfile with filename "{bamfile}" '
                                     f'can be found. Are you sure you have run the pipeline mapping step already for '
                                     f'this sample?', width=90, subsequent_indent=' ' * 11)
                logger.info(fill)
                sys.exit(1)
        elif args.blast:
            blastx_outputfile = f'{basename}.blastx'
            logger.debug(f'blastx outputfile is: {blastx_outputfile}')
            if not utils.file_exists_and_not_empty(blastx_outputfile):
                fill = textwrap.fill(f'{"[ERROR]:":10} The parameter "--start_from {args.start_from}" has been provided '
                                     f' and BLASTx/DIAMOND is used for read mapping, but no existing *.blastx output '
                                     f'file with filename "{blastx_outputfile}" can be found. Are you sure you have '
                                     f'run the pipeline mapping step already for this sample?', width=90,
                                     subsequent_indent=' ' * 11)
                logger.info(fill)
                sys.exit(1)
    else:
        if args.bwa:  # map reads to nucleotide targets with BWA
            if args.unpaired:
                # Note that unpaired_readfile is a single path to the file:
                bwa(unpaired_readfile, targetfile, basename, cpu=cpu, unpaired=True, logger=logger)
            # Note that readfiles is a list of one (single-end) or two (paired-end) paths to read files:
            bamfile = bwa(readfiles, targetfile, basename, cpu=cpu, logger=logger)
            if not bamfile:
                logger.error(f'{"[ERROR]:":10} Something went wrong with the BWA step, exiting. Check the '
                             f'hybpiper_assemble.log file for sample {basename}!')
                return
            logger.debug(f'bamfile is: {bamfile}')

        elif args.blast:  # map reads to protein targets with BLASTx
            if args.unpaired:
                blastx(unpaired_readfile, targetfile, args.evalue, basename, cpu=cpu,
                       max_target_seqs=args.max_target_seqs, unpaired=True, logger=logger, diamond=args.diamond,
                       diamond_sensitivity=args.diamond_sensitivity)

            blastx_outputfile = blastx(readfiles, targetfile, args.evalue, basename, cpu=cpu,
                                       max_target_seqs=args.max_target_seqs, logger=logger, diamond=args.diamond,
                                       diamond_sensitivity=args.diamond_sensitivity)

            if not blastx_outputfile:
                logger.error(f'{"[ERROR]:":10} Something went wrong with the Blastx step, exiting. Check the '
                             f'hybpiper_assemble.log file for sample {basename}!')
                return
        else:
            sys.exit(f'Can not determine whether BWA or BLASTx option is supplied, exiting...')

    ####################################################################################################################
    # Distribute reads to gene directories from either BLASTx or BWA mapping
    ####################################################################################################################

    if args.start_from not in ['map_reads', 'distribute_reads']:
        logger.info(f'{"[INFO]:":10} Parameter "--start_from {args.start_from}" supplied, skipping read and target '
                    f'distribution step!')
        pre_existing_fastas = glob.glob('./*/*_interleaved.fasta') + glob.glob('./*/*_unpaired.fasta')
        if len(pre_existing_fastas) == 0:
            fill = textwrap.fill(f'{"[ERROR]:":10} The parameter "--start_from {args.start_from}" has been provided but '
                                 f'no distributed reads (*_interleaved.fasta and/or *_unpaired.fasta) can be found '
                                 f'for any gene. Are you sure you have run the pipeline read distribution step '
                                 f'already for this sample?', width=90,
                                 subsequent_indent=' ' * 11)
            logger.info(fill)
            sys.exit(1)
    else:
        pre_existing_fastas = glob.glob('./*/*_interleaved.fasta') + glob.glob('./*/*_unpaired.fasta')
        for fasta in pre_existing_fastas:
            os.remove(fasta)

        if args.bwa:
            distribute_bwa(bamfile, readfiles, targetfile, target, unpaired_readfile, args.exclude,
                           merged=args.merged, hi_mem=args.distribute_hi_mem, logger=logger)
        else:  # distribute BLASTx results
            distribute_blastx(blastx_outputfile, readfiles, targetfile, target, unpaired_readfile, args.exclude,
                              merged=args.merged, hi_mem=args.distribute_hi_mem, logger=logger)

    # Note that HybPiper expects either paired-end readfiles (parameter --readfiles) and an optional file of unpaired
    #  reads (parameter --unpaired), or a single file of unpaired reads (parameter --readfiles). For each scenario,
    #  the unpaired readfile is written to an *_unpaired.fasta file.
    if len(readfiles) == 2:
        genes = [x for x in os.listdir('.') if os.path.isfile(os.path.join(x, x + '_interleaved.fasta'))]
    else:
        genes = [x for x in os.listdir('.') if os.path.isfile(os.path.join(x, x + '_unpaired.fasta'))]
    if len(genes) == 0:
        logger.error('ERROR: No genes with reads, exiting!')
        return

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

    ####################################################################################################################
    # Run Exonerate on the assembled SPAdes contigs, and intronerate() if flag --run_intronerate is supplied:
    ####################################################################################################################

    genes = [x.rstrip() for x in open('exonerate_genelist.txt').readlines()]
    # Remove any pre-existing directories:
    for g in genes:
        if os.path.isdir(os.path.join(g, basename)):
            shutil.rmtree(os.path.join(g, basename))
    if os.path.isfile('genes_with_seqs.txt'):
        os.remove('genes_with_seqs.txt')
    if len(genes) == 0:
        logger.error(f'{"[ERROR]:":10} No genes recovered for {basename}!')
        return 1

    exonerate_multiprocessing(genes,
                              basename,
                              thresh=args.thresh,
                              paralog_warning_min_length_percentage=args.paralog_min_length_percentage,
                              depth_multiplier=args.depth_multiplier,
                              no_stitched_contig=args.no_stitched_contig,
                              bbmap_memory=args.bbmap_memory,
                              bbmap_subfilter=args.bbmap_subfilter,
                              chimeric_stitched_contig_edit_distance=args.chimeric_stitched_contig_edit_distance,
                              chimeric_stitched_contig_discordant_reads_cutoff=
                              args.chimeric_stitched_contig_discordant_reads_cutoff,
                              bbmap_threads=args.bbmap_threads,
                              pool_threads=cpu,
                              logger=logger,
                              intronerate=args.intronerate,
                              no_padding_supercontigs=args.no_padding_supercontigs,
                              keep_intermediate_files=args.keep_intermediate_files,
                              exonerate_contigs_timeout=args.timeout_exonerate_contigs,
                              verbose_logging=args.verbose_logging)

    ####################################################################################################################
    # Collate all stitched contig and putative chimera read reports
    ####################################################################################################################
    logger.info(f'\n{"[INFO]:":10} Generated sequences from {len(open("genes_with_seqs.txt").readlines())} genes!')

    # Stitched contigs:
    collate_stitched_contig_reports = [x for x in glob.glob(f'*/{basename}/genes_with_stitched_contig.csv')]
    with open(f'{basename}_genes_with_stitched_contig.csv', 'w') as genes_with_stitched_contig_handle:
        for report_file in collate_stitched_contig_reports:
            with open(report_file, 'r') as report_handle:
                lines = report_handle.readlines()
                genes_with_stitched_contig_handle.write('\n'.join(lines))

    # Putative chimeras:
    collate_putative_chimeras_reports = [x for x in glob.glob(f'*/{basename}/putative_chimeric_stitched_contig.csv')]
    with open(f'{basename}_genes_derived_from_putative_chimeric_stitched_contig.csv',
              'w') as genes_with_chimeras_handle:
        for report_file in collate_putative_chimeras_reports:
            with open(report_file, 'r') as report_handle:
                lines = report_handle.readlines()
                genes_with_chimeras_handle.write('\n'.join(lines))

    ####################################################################################################################
    # Report paralog warnings and write paralog warning files
    ####################################################################################################################

    # Collate report for long paralogs, and write warning to screen:
    paralog_warnings_long = [x for x in glob.glob(f'*/{basename}/paralog_warning_long.txt')]
    with open(f'{basename}_genes_with_long_paralog_warnings.txt', 'w') as long_paralogs_handle:
        for warning_file in paralog_warnings_long:
            with open(warning_file, 'r') as warning_handle:
                report_line = warning_handle.readline().split()[0]  # only recover gene name
                long_paralogs_handle.write(f'{report_line}\n')
    logger.info(f'{"[WARNING]:":10} Potential long paralogs detected for {len(paralog_warnings_long)} genes!')

    # Collate report for paralogs via SPAdes contig depth, and write warning to screen:
    paralog_warnings_short = [x for x in glob.glob(f'*/{basename}/paralog_warning_by_contig_depth.txt')]
    paralog_warnings_short_true = 0
    with open(f'{basename}_genes_with_paralog_warnings_by_contig_depth.csv', 'w') as depth_paralogs_handle:
        for warning_file in paralog_warnings_short:
            with open(warning_file, 'r') as warning_handle:
                report_line = warning_handle.readline()
                if report_line.split()[-1] == 'True':
                    paralog_warnings_short_true += 1
                depth_paralogs_handle.write(report_line)
    logger.info(f'{"[WARNING]:":10} Potential paralogs detected via contig depth for'
                f' {paralog_warnings_short_true} genes!')

    logger.info(f'\nFinished running "hybpiper assemble" for sample {basename}!\n')


def hybpiper_stats_main(args):
    """
    Calls the function main() from module hybpiper_stats

    :param args: argparse namespace with subparser options for function get_seq_lengths_main()
    :return: None: no return value specified; default is None
    """

    hybpiper_stats.main(args)


def retrieve_sequences_main(args):
    """
    Calls the function main() from module retrieve_sequences

    :param args: argparse namespace with subparser options for function retrieve_sequences_main()
    :return: None: no return value specified; default is None
    """

    retrieve_sequences.main(args)


def paralog_retriever_main(args):
    """
    Calls the function main() from module paralog_retriever

    :param args: argparse namespace with subparser options for function paralog_retriever_main()
    :return: None: no return value specified; default is None
    """

    paralog_retriever.main(args)


def gene_recovery_heatmap_main(args):
    """
    Calls the function main() from module gene_recovery_heatmap

    :param args: argparse namespace with subparser options for function gene_recovery_heatmap_main()
    :return: None: no return value specified; default is None
    """

    gene_recovery_heatmap.main(args)


def check_dependencies(args):
    """
    # Calls the function check_dependencies() from module utils

    :param args: argparse namespace with subparser options for function check_dependencies()
    :return: None: no return value specified; default is None
    """

    utils.check_dependencies(logger=args.logger)


def check_targetfile_standalone(args):
    """
    Performs targetfile checks. Does not translate a DNA file; low-complexity checks are performed on the target file
    as provided.

    :param args: argparse namespace with subparser options for function check_targetfile()
    :return: None: no return value specified; default is None
    """

    # Set target file type and path:
    if args.targetfile_dna:
        targetfile = args.targetfile_dna
        targetfile_type = 'DNA'
        translate_target_file = True
    elif args.targetfile_aa:
        targetfile = args.targetfile_aa
        targetfile_type = 'protein'
        translate_target_file = False

    # Check targetfile header and duplicate gene names:
    check_target_file_headers_and_duplicate_names(targetfile, logger=args.logger)

    # Check that seqs in target file can be translated from the first codon position in the forwards frame:
    check_target_file_stop_codons_and_multiple_of_three(targetfile, translate_target_file=translate_target_file,
                                                        logger=None)
    low_complexity_sequences = None

    if targetfile_type == 'DNA':
        fill = textwrap.fill(f'{"[INFO]:":10} The target file {targetfile} has been specified as containing DNA '
                             f'sequences. These DNA sequences will be checked for low-complexity regions. NOTE: if you '
                             f'run "hybpiper assemble" without providing the flag "--bwa" the DNA target file will be '
                             f'translated and the check will be run again on the translated protein sequences; the '
                             f'sequences flagged as having low-complexity regions can sometimes differ between a DNA '
                             f'and a corresponding translated protein target file.', width=90,
                             subsequent_indent=" " * 11)
        print(fill)

        low_complexity_sequences = utils.low_complexity_check(targetfile,
                                                              targetfile_type,
                                                              translate_target_file=False,
                                                              window_size=args.sliding_window_size,
                                                              entropy_value=args.complexity_minimum_threshold,
                                                              logger=args.logger)
    elif targetfile_type == 'protein':
        fill = textwrap.fill(f'{"[INFO]:":10} The target file {targetfile} has been specified as containing protein '
                             f'sequences. These protein sequences will be checked for low-complexity regions', width=90,
                             subsequent_indent=" " * 11)
        print(fill)

        low_complexity_sequences = utils.low_complexity_check(targetfile,
                                                              targetfile_type,
                                                              translate_target_file=False,
                                                              window_size=args.sliding_window_size,
                                                              entropy_value=args.complexity_minimum_threshold,
                                                              logger=args.logger)

    if low_complexity_sequences:
        fill_1 = textwrap.fill(f'{"[WARNING]:":10} The target file provided ({os.path.basename(targetfile)}) contains '
                               f'sequences with low-complexity regions. The sequence names are printed below. These '
                               f'sequences can cause problems when running HybPiper, '
                               f'see https://github.com/mossmatters/HybPiper/wiki/Troubleshooting,-common-issues,'
                               f'-and-recommendations. We recommend one of the following approaches:', width=90,
                               subsequent_indent=" " * 11)

        fill_2 = textwrap.fill(f'1) Remove these sequence from your target file, ensuring that your file still '
                               f'contains other representative sequences for the corresponding genes.', width=90,
                               initial_indent=" " * 11, subsequent_indent=" " * 14)

        fill_3 = textwrap.fill(f'2) Start the run using the parameter "--timeout_assemble" (e.g. "--timeout_assemble '
                               f'200"). See '
                               f'https://github.com/mossmatters/HybPiper/wiki/Full-pipeline-parameters#10-hybpiper'
                               f'-assemble for details.',
                               width=90, initial_indent=" " * 11, subsequent_indent=" " * 14, break_on_hyphens=False)

        print(f'{fill_1}\n\n{fill_2}\n\n{fill_3}\n')
        print(f'\n{" " * 10} Sequences with low complexity regions are:\n')

        for sequence in low_complexity_sequences:
            print(f'{" " * 10} {sequence}')

    else:
        print(f'{"[INFO]:":10} No sequences with low-complexity regions found.')


def parse_arguments():
    """
    Creates main parser and add subparsers. Parses command line arguments

    :return argparse.Namespace arguments: arguments for the given command/subcommand
    """

    parser = argparse.ArgumentParser(prog='hybpiper', description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='To view parameters and help for a subcommand, use e.g. "assemble '
                                            '--help"')
    group_1 = parser.add_mutually_exclusive_group(required=False)
    group_1.add_argument('--version', '-v',
                         dest='version',
                         action='version',
                         version='%(prog)s 2.0.1rc build 8',
                         help='Print the HybPiper version number.')

    # Add subparsers:
    subparsers = parser.add_subparsers(title='Subcommands for HybPiper', description='Valid subcommands:')
    parser_assemble = hybpiper_subparsers.add_assemble_parser(subparsers)
    parser_stats = hybpiper_subparsers.add_stats_parser(subparsers)
    parser_retrieve_sequences = hybpiper_subparsers.add_retrieve_sequences_parser(subparsers)
    parser_paralog_retriever = hybpiper_subparsers.add_paralog_retriever_parser(subparsers)
    parser_gene_recovery_heatmap = hybpiper_subparsers.add_gene_recovery_heatmap_parser(subparsers)
    parser_check_dependencies = hybpiper_subparsers.add_check_dependencies_parser(subparsers)
    parser_check_targetfile = hybpiper_subparsers.add_check_targetfile_parser(subparsers)

    # Set functions for subparsers:
    parser_assemble.set_defaults(func=assemble)
    parser_stats.set_defaults(func=hybpiper_stats_main)
    parser_retrieve_sequences.set_defaults(func=retrieve_sequences_main)
    parser_paralog_retriever.set_defaults(func=paralog_retriever_main)
    parser_gene_recovery_heatmap.set_defaults(func=gene_recovery_heatmap_main)
    parser_check_dependencies.set_defaults(func=check_dependencies)
    parser_check_targetfile.set_defaults(func=check_targetfile_standalone)

    # Parse and return all arguments:
    arguments = parser.parse_args()

    # Get the run directory containing the assemble.py module:
    run_dir = os.path.realpath(os.path.split(sys.argv[0])[0])

    return arguments


def main():

    if len(sys.argv) == 1:
        print(__doc__)
        sys.exit(1)

    # Parse arguments for the command/subcommand used:
    args = parse_arguments()

    # Run the function associated with the subcommand:
    args.func(args)


########################################################################################################################
# Run the script
#######################################################################################################################
if __name__ == '__main__':
    main()
    # import pstats
    # import cProfile
    # profiler = cProfile.Profile()
    # profiler.enable()
    # main()
    # profiler.disable()
    # stats = pstats.Stats(profiler).sort_stats('cumtime')
    # stats.print_stats()

################################################## END OF SCRIPT #######################################################
