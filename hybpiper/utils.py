#!/usr/bin/env python

"""
This module contains some general functions and classes.
"""
import copy
import re
from textwrap import TextWrapper
import os
import collections
import scipy
import textwrap
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Seq import Seq
import subprocess
import datetime
import logging
import sys
import io
import pstats
import cProfile
import platform
import resource
import argparse
import progressbar
from collections import defaultdict
import glob
import tarfile
import shutil

from hybpiper.version import __version__


def log_or_print(string, logger=None, logger_level='info'):
    """
    If a logger is provided, print to string using logger. If not, use print()

    :param logging.Logger logger: a logger object
    :param str string: a string to print
    :param str logger_level: level to use for logger e.g. info/debug/error
    :return:
    """

    if logger:
        if logger_level == 'info':
            logger.info(string)
        elif logger_level == 'debug':
            logger.debug(string)
        elif logger_level == 'error':
            logger.error(string)
        elif logger_level == 'warning':
            logger.warning(string)
    else:
        print(string)


def file_exists_and_not_empty(file_name):
    """
    Check if file exists and is not empty by confirming that its size is not 0 bytes. Returns a boolean.

    :param str file_name: path to filename to check
    :return: bool
    """

    return os.path.isfile(file_name) and not os.path.getsize(file_name) == 0


def fill_forward_slash(text, width=70, **kwargs):
    """
    Fill a single paragraph of text, returning a new string.

    Reformat the single paragraph in 'text' to fit in lines of no more
    than 'width' columns, and return a new string containing the entire
    wrapped paragraph.  As with wrap(), tabs are expanded and other
    whitespace characters converted to space.  See TextWrapper class for
    available keyword args to customize wrapping behaviour.

    This function uses the subclass TextWrapperForwardSlash.
    """
    w = TextWrapperForwardSlash(width=width, **kwargs)
    return w.fill(text)


class TextWrapperForwardSlash(TextWrapper):
    """
    Subclasses textwrap.TextWrapper and alters regex so that break_on_hyphens corresponds to forward slashes rather
    than hyphens. Used for wrapping long path strings.

    Change: letter = r'[^\d\W]' -> letter = r'[\w-]'
    """

    _whitespace = '\t\n\x0b\x0c\r '
    word_punct = r'[\w!"\'&.,?]'
    letter = r'[\w-]'
    whitespace = r'[%s]' % re.escape(_whitespace)
    nowhitespace = '[^' + whitespace[1:]

    wordsep_re = re.compile(r'''
        ( # any whitespace
          %(ws)s+
        | # em-dash between words
          (?<=%(wp)s) -{2,} (?=\w)
        | # word, possibly hyphenated
          %(nws)s+? (?:
            # hyphenated word
              (?: (?<=%(lt)s{2}/) | (?<=%(lt)s/%(lt)s/))
              (?= %(lt)s /? %(lt)s)
            | # end of word
              (?=%(ws)s|\Z)
            | # em-dash
              (?<=%(wp)s) (?=-{2,}\w)
            )
        )''' % {'wp': word_punct,
                'lt': letter,
                'ws': whitespace,
                'nws': nowhitespace},
        re.VERBOSE)
    del word_punct, letter, nowhitespace

    def __init__(self,
                 width=70,
                 initial_indent="",
                 subsequent_indent="",
                 expand_tabs=True,
                 replace_whitespace=True,
                 fix_sentence_endings=False,
                 break_long_words=True,
                 drop_whitespace=True,
                 break_on_forward_slash=True,
                 tabsize=8,
                 *,
                 max_lines=None,
                 placeholder=' [...]'):
        super().__init__(width=width,
                         initial_indent=initial_indent,
                         subsequent_indent=subsequent_indent,
                         expand_tabs=expand_tabs,
                         replace_whitespace=replace_whitespace,
                         fix_sentence_endings=fix_sentence_endings,
                         break_long_words=break_long_words,
                         drop_whitespace=drop_whitespace,
                         break_on_hyphens=break_on_forward_slash,
                         tabsize=tabsize,
                         max_lines=max_lines,
                         placeholder=placeholder)


def py_which(executable_name, mode=os.F_OK | os.X_OK, path=None):
    """
    Given an executable_name, mode, and a PATH string, return the path which conforms to the given mode on the PATH,
    or None if there is no such file. `mode` defaults to os.F_OK | os.X_OK. `path` defaults to the result
    of os.environ.get("PATH"), or can be overridden with a custom search path.

    :param str executable_name: executable name to search for
    :param int mode: bitwise OR of integers from os.F_OK (existence of path) and os.X_OK (executable) access checks
    :param str path: None, or contents of $PATH variable (as recovered in function body)
    :return None or path of executable
    """

    # Check that a given file can be accessed with the correct mode. Additionally, check that `file` is not a
    # directory, as on Windows directories pass the os.access check.
    def _access_check(fn, mode):
        return (os.path.exists(fn) and os.access(fn, mode)
                and not os.path.isdir(fn))

    # If we're given a path with a directory part, look it up directly rather than referring to PATH directories.
    # This includes checking relative to the current directory, e.g. ./script
    if os.path.dirname(executable_name):
        if _access_check(executable_name, mode):
            return executable_name
        return None

    if path is None:
        path = os.environ.get('PATH', os.defpath)
    if not path:
        return None

    path = path.split(os.pathsep)
    checked_directories = set()

    for directory in path:
        if directory not in checked_directories:  # only check each directory once
            checked_directories.add(directory)
            name = os.path.join(directory, executable_name)
            if _access_check(name, mode):
                return name

    return None


def shannon_entropy(fasta_seq):
    """
    Calculate Shannon entropy for a given string.

    :param fasta_seq: a string of characters (e.g. DNA or amino-acid sequence).
    :return: shannon entropy value (float) for the given string.
    """

    characters = collections.Counter([tmp_character for tmp_character in fasta_seq])
    # e.g. Counter({'A': 4, 'T': 4, 'G': 1,'C': 1})
    dist = [x / sum(characters.values()) for x in characters.values()]  # e.g. [0.4, 0.4, 0.1, 0.1]

    # use scipy to calculate entropy
    entropy_value = scipy.stats.entropy(dist, base=2)

    return entropy_value


def low_complexity_check(targetfile, targetfile_type, translate_target_file, window_size=None, entropy_value=None,
                         logger=None):
    """
    For each sequence in a fasta file, recover the sequence string of a sliding window and calculate the Shannon
    entropy value for it. If it is below a threshold, flag the entire sequence as low entropy and add it to a set.
    Entropy will be calculated from protein sequences if a nucleotide target file has been translated for BLASTx or
    DIAMOND. Returns the set of sequence names.

    :param str targetfile: path to the targetfile
    :param str targetfile_type: string describing target file sequence type i.e 'DNA' or 'protein'
    :param bool translate_target_file: True is a nucleotide target file has been translated for BLASTx orDIAMOND
    :param int or None window_size: number of characters to include in the sliding window
    :param float or None entropy_value: entropy threshold within sliding window
    :param logging.Logger logger: a logger object
    :return set low_entropy_seqs, int window_size, float entropy_value: a set of sequences that contain low entropy
    substrings, an integer corresponding to sliding window size, a float corresponding to the entropy value threshold
    """

    # Set widget format for progressbar:
    widgets = [' ' * 11,
               progressbar.Timer(),
               progressbar.Bar(),
               progressbar.ETA()]

    if targetfile_type == 'DNA' and not translate_target_file:
        entropy_value = entropy_value if entropy_value else 1.5
        window_size = window_size if window_size else 100
    elif targetfile_type == 'protein' or (targetfile_type == 'DNA' and translate_target_file):
        entropy_value = entropy_value if entropy_value else 3.0
        window_size = window_size if window_size else 50

    fill_1 = textwrap.fill(f'{"[INFO]:":10} Checking the target file for sequences with low-complexity regions...',
                           width=90, subsequent_indent=" " * 11)
    fill_2 = textwrap.fill(f'{"[INFO]:":10} Type of sequences: {targetfile_type}', width=90,
                           subsequent_indent=" " * 11)
    fill_3 = textwrap.fill(f'{"[INFO]:":10} Sliding window size: {window_size}', width=90,
                           subsequent_indent=" " * 11)
    fill_4 = textwrap.fill(f'{"[INFO]:":10} Minimum complexity threshold: {entropy_value}', width=90,
                           subsequent_indent=" " * 11)

    log_or_print(fill_1, logger=logger)
    log_or_print(fill_2, logger=logger)
    log_or_print(fill_3, logger=logger)
    log_or_print(fill_4, logger=logger)

    total_seqs = [seq for seq in SeqIO.parse(targetfile, "fasta")]

    low_entropy_seqs = set()
    for sequence in progressbar.progressbar(total_seqs, max_value=len(total_seqs), min_poll_interval=30,
                                            widgets=widgets):

        for i in range(0, len(sequence.seq) - (window_size - 1)):
            window_seq = str(sequence.seq[i:i + window_size])
            window_shannon_entropy = shannon_entropy(window_seq)
            if window_shannon_entropy <= entropy_value:
                low_entropy_seqs.add(sequence.name)
                break

    if len(low_entropy_seqs) == 0:
        log_or_print(f'{"[INFO]:":10} No sequences with low-complexity regions found.', logger=logger)

    return low_entropy_seqs, window_size, entropy_value


def pad_seq(sequence):
    """
    Pads a sequence Seq object to a multiple of 3 with 'N'.

    :param Bio.SeqRecord.SeqRecord sequence: sequence to pad
    :return: Bio.SeqRecord.SeqRecord sequence padded with Ns if required, bool True if sequence needed padding
    """

    remainder = len(sequence.seq) % 3
    if remainder == 0:
        return sequence, False
    else:
        # logger.info(f'{"[WARN!]:":10} The targetfile nucleotide sequence {sequence.id} is not a multiple of 3!')
        sequence.seq = sequence.seq + Seq('N' * (3 - remainder))
        return sequence, True


def check_exonerate_version():
    """
    Returns the Exonerate version e.g. 2.4.0

    :return str version: the Exonerate version number
    """

    exonerate_help = subprocess.run(['exonerate', '-h'], capture_output=True)
    version = re.search(r'\bexonerate version [0-9][.][0-9].[0-9]\b', str(exonerate_help)).group(0)
    version_number = version.split()[-1]
    return version_number


def check_dependencies(logger=None):
    """
    Checks for the presence of executables and Python packages. Returns a boolean.

    :param None, logging.Logger logger: a logger object
    return: bool everything_is_awesome: True if all dependencies are found and are executable, else False
    """

    executables = ['blastx',
                   'exonerate',
                   'parallel',
                   'makeblastdb',
                   'spades.py',
                   'bwa',
                   'samtools',
                   'bbmap.sh',
                   'bbmerge.sh',
                   'diamond',
                   'mafft']

    log_or_print(f'{"[INFO]:":10} Checking for external dependencies:\n', logger=logger)

    everything_is_awesome = True
    for exe in executables:
        exe_loc = py_which(exe)
        if exe_loc:
            log_or_print(f'{exe:20} found at {exe_loc}', logger=logger)
            if exe == 'exonerate':
                version = check_exonerate_version()
                if version not in ['2.4.0']:
                    everything_is_awesome = False

                    log_or_print(f'\n{exe} version 2.4.0 is required, but your version is {version}. Please '
                                 f'update your Exonerate version!\n', logger=logger)
        else:
            log_or_print(f'{exe:20} not found in your $PATH!', logger=logger)
            everything_is_awesome = False

    log_or_print('', logger=logger)

    return everything_is_awesome


def make_basename(readfiles, prefix=None, output_folder=None):
    """
    Unless prefix is set, generate a directory name based on the readfiles name, using everything up to the first
    underscore. If prefix is set, generate the directory name "prefix". In both cases, if --hybpiper_output is set,
    generate the full path to the sample folder using the specified parent output directory.

    https://docs.python.org/3/library/os.path.html#os.path.split:

    Split the pathname path into a pair, (head, tail) where tail is the last pathname component and head is
    everything leading up to that. The tail part will never contain a slash; if path ends in a slash, tail will be
    empty. If there is no slash in path, head will be empty. If path is empty, both head and tail are empty. Trailing
    slashes are stripped from head unless it is the root (one or more slashes only). In all cases, join(head, tail)
    returns a path to the same location as path (but the strings may differ). Also see the functions dirname() and
    basename().

    :param list readfiles: one or more read files used as input to the pipeline
    :param str prefix: directory name for sample pipeline output
    :param str/None output_folder: name of output folder if supplied, else None
    :return str parent directory, directory name
    """

    parent_dir = os.path.abspath(output_folder) if output_folder else os.getcwd()

    if prefix:
        full_sample_directory = os.path.join(parent_dir, prefix)

    else:  # If --prefix is not set on cmd line, get sample name from read file:
        prefix = os.path.split(readfiles[0])[1].split('_')[0]
        full_sample_directory = os.path.join(parent_dir, prefix)

    return (full_sample_directory,
            parent_dir,
            prefix)


def worker_configurer(gene_name):
    """
    Configures logging to file and screen for the worker processes

    :param str gene_name: name of the gene being processing by the worker process
    :return: None
    """

    # Get date and time string for log filename:
    date_and_time = datetime.datetime.now().strftime("%Y-%m-%d-%H_%M_%S")

    # Log to file:
    file_handler = logging.FileHandler(f'{gene_name}/{gene_name}_{date_and_time}.log', mode='w')
    file_handler.setLevel(logging.DEBUG)
    file_format = logging.Formatter('%(asctime)s - %(filename)s - %(name)s - %(funcName)s - %(levelname)s - %('
                                    'message)s')
    file_handler.setFormatter(file_format)

    # Log to Terminal (stdout):
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(logging.INFO)
    console_format = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_format)

    # Setup logger:
    logger_object = logging.getLogger(gene_name)

    # Add handlers to the logger
    logger_object.addHandler(console_handler)
    logger_object.addHandler(file_handler)


def cprofile_to_csv(profile_binary_file):
    """
    Takes a cProfile.Profile object, converts it to human-readable format, and returns the data in *.csv format for
    writing to file.
    From: https://gist.github.com/ralfstx/a173a7e4c37afa105a66f371a09aa83e
    :param cProfile.Profile profile_binary_file: a cProfile.Profile object
    :return str: human-readable data from the cProfile run, in *.csv format
    """

    out_stream = io.StringIO()
    pstats.Stats(profile_binary_file, stream=out_stream).sort_stats('cumtime').print_stats()
    result = out_stream.getvalue()
    result = 'ncalls' + result.split('ncalls')[-1]  # chop off header lines
    lines = [','.join(line.rstrip().split(None, 5)) for line in result.split('\n')]

    return '\n'.join(lines)


def write_fix_targetfile_controlfile(targetfile_type,
                                     translate_target_file,
                                     no_terminal_stop_codons,
                                     low_complexity_sequences,
                                     sliding_window_size,
                                     complexity_minimum_threshold):
    """
    Write a control file with settings and seqs to remove, as output from the command `hybpiper check_targetfile`.

    :param str targetfile_type:
    :param bool translate_target_file:
    :param bool no_terminal_stop_codons:
    :param lists / None low_complexity_sequences:
    :param int sliding_window_size:
    :param float complexity_minimum_threshold:
    """

    date_and_time = datetime.datetime.now().strftime("%Y-%m-%d-%H_%M_%S")

    low_complexity_sequences_to_remove = '\t'.join(low_complexity_sequences) if low_complexity_sequences else None

    with open(f'fix_targetfile_{date_and_time}.ctl', 'w') as targetfile_control_handle:
        targetfile_control_handle.write(f'TARGETFILE_TYPE\t{targetfile_type}\n')
        targetfile_control_handle.write(f'TRANSLATE_TARGET_FILE\t{translate_target_file}\n')
        targetfile_control_handle.write(f'NO_TERMINAL_STOP_CODONS\t{no_terminal_stop_codons}\n')
        targetfile_control_handle.write(f'SLIDING_WINDOW_SIZE\t{sliding_window_size}\n')
        targetfile_control_handle.write(f'COMPLEXITY_MINIMUM_THRESHOLD\t{complexity_minimum_threshold}\n')
        targetfile_control_handle.write(f'ALLOW_GENE_REMOVAL\tFalse\n')
        targetfile_control_handle.write(f'LOW_COMPLEXITY_SEQUENCES\t{low_complexity_sequences_to_remove}\n')

    print(f'\n{"[INFO]:":10} The control file required for command "hybpiper fix_targetfile" has been written to\n'
          f'           "fix_targetfile_{date_and_time}.ctl".')


def restricted_float(input_float):
    """
    Checks that a provided value is a float within the range 0.0 to 1.0

    From: https://stackoverflow.com/questions/12116685/\
    how-can-i-require-my-python-scripts-argument-to-be-a-float-in
    -a-range-using-arg
    """

    try:
        input_float = float(input_float)
    except ValueError:
        raise argparse.ArgumentTypeError(f'{input_float} not a floating-point literal')

    if input_float < 0.0 or input_float > 1.0:
        raise argparse.ArgumentTypeError(f'{input_float} not in range [0.0, 1.0]')

    return input_float


def restricted_int_word_size(input_int,
                             start_value=4,
                             end_value=1000):
    """
    Checks that a provided value is an integer within the range start_value to end_value

    :param input_int:
    :param int start_value:
    :param int end_value:
    :return:
    """

    try:
        input_int = int(input_int)
    except ValueError:
        raise argparse.ArgumentTypeError(f'{input_int} not an integer!')

    if input_int < start_value or input_int > end_value:
        raise argparse.ArgumentTypeError(f'{input_int} not in range [{start_value}, {end_value}]')

    return input_int


def restricted_int_perc_identity(input_int,
                                 start_value=0,
                                 end_value=100):
    """
    Checks that a provided value is an integer within the range start_value to end_value

    :param input_int:
    :param int start_value:
    :param int end_value:
    :return:
    """

    try:
        input_int = int(input_int)
    except ValueError:
        raise argparse.ArgumentTypeError(f'{input_int} not an integer!')

    if input_int < start_value or input_int > end_value:
        raise argparse.ArgumentTypeError(f'{input_int} not in range [{start_value}, {end_value}]')

    return input_int


def createfolder(directory):
    """
    Attempts to create a directory named after the name provided, and provides an error message on failure
    """
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
        return directory
    except OSError:
        print(f'Error: Creating directory: {directory}')
        raise

def check_macos_version(logger=None):
    """
    Due to this issue:

        https://stackoverflow.com/questions/65290242/pythons-platform-mac-ver-reports-incorrect-macos-version

    ...use the macOS command `sw_vers` to get the correct macOS version.

    :param None, logging.Logger logger: a logger object
    """

    try:
        result = subprocess.run('sw_vers', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True, check=True)
        logger.debug(f'sw_vers check_returncode() is: {result.check_returncode()}')
        logger.debug(f'sw_vers stdout is: {result.stdout}')
        logger.debug(f'sw_vers stderr is: {result.stderr}')

    except subprocess.CalledProcessError as exc:
        logger.error(f'sw_vers FAILED. Output is: {exc}')
        logger.error(f'sw_vers stdout is: {exc.stdout}')
        logger.error(f'sw_vers stderr is: {exc.stderr}')


def get_platform_info(logger=None):
    """
    Log the platform version for debugging

    :param None, logging.Logger logger: a logger object
    """

    logger.debug(f'uname:     {platform.uname()}')
    logger.debug(f'system:    {platform.system()}')
    logger.debug(f'node:      {platform.node()}')
    logger.debug(f'release:   {platform.release()}')
    logger.debug(f'version:   {platform.version()}')
    logger.debug(f'machine:   {platform.machine()}')
    logger.debug(f'processor: {platform.processor()}')


def get_ulimit_info(logger=None):
    """
    Log ulimit details for debugging.

    Uses: https://github.com/python/cpython/blob/main/Modules/resource.c

    :param None, logging.Logger logger: a logger object
    """

    for name, desc in [
        ('RLIMIT_CORE', 'core file size'),
        ('RLIMIT_CPU', 'CPU time'),
        ('RLIMIT_FSIZE', 'file size'),
        ('RLIMIT_DATA', 'heap size'),
        ('RLIMIT_STACK', 'stack size'),
        ('RLIMIT_RSS', 'resident set size'),
        ('RLIMIT_NPROC', 'number of processes'),
        ('RLIMIT_NOFILE', 'number of open files'),
        ('RLIMIT_MEMLOCK', 'lockable memory address'),
    ]:
        try:
            limit_num = getattr(resource, name)
            soft, hard = resource.getrlimit(limit_num)
            logger.debug(f'Maximum {desc:25} ({name:15}) : {soft:20} {hard:20}')
        except ValueError:
            logger.info(f'Specified resource {name} not found!')


def check_target_file_headers_and_duplicate_names(targetfile, logger=None):
    """
    - Checks target-file fasta header formatting ("taxon*-unique_gene_ID").
    - Checks for duplicate gene names in the targetfile.
    - Reports the number of unique genes (each can have multiple representatives) in the targetfile.

    :param str targetfile: path to the targetfile
    :param logging.Logger logger: a logger object
    :return: incorrectly_formatted_fasta_headers, duplicated_gene_names, gene_lists
    """

    gene_lists = defaultdict(list)
    with open(targetfile, 'r') as target_file_handle:
        seqs = list(SeqIO.parse(target_file_handle, 'fasta'))
        incorrectly_formatted_fasta_headers = set()
        check_for_duplicate_genes_dict = {}
        for seq in seqs:
            if seq.name in check_for_duplicate_genes_dict:
                check_for_duplicate_genes_dict[seq.name] += 1
            else:
                check_for_duplicate_genes_dict[seq.name] = 1
            if not re.match('.+-[^-]+$', seq.name):
                incorrectly_formatted_fasta_headers.add(seq.name)
            if re.search('\"|\'', seq.name):
                incorrectly_formatted_fasta_headers.add(seq.name)
            gene_id = re.split('-', seq.name)[-1]
            gene_lists[gene_id].append(seq)

    if len(incorrectly_formatted_fasta_headers) == 0:
        log_or_print(f'{"[INFO]:":10} The target file FASTA header formatting looks good!', logger=logger)

    # Check for duplicated gene names:
    duplicated_gene_names = []
    for gene, gene_count in check_for_duplicate_genes_dict.items():
        if gene_count > 1:
            duplicated_gene_names.append(gene)

    # Report the number of unique genes represented in the target file:
    # log_or_print(f'{"[INFO]:":10} The target file contains at least one sequence for {len(gene_lists)} '
    #              f'unique genes.', logger=logger)

    return incorrectly_formatted_fasta_headers, duplicated_gene_names, gene_lists


def setup_logger(name, log_file, console_level=logging.INFO, file_level=logging.DEBUG,
                 logger_object_level=logging.DEBUG):
    """
    Function to create a logger instance.

    By default, logs level DEBUG and above to file.
    By default, logs level INFO and above to stdout and file.

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
    console_handler = logging.StreamHandler(sys.stdout)
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


def check_target_file_stop_codons_and_multiple_of_three(targetfile,
                                                        no_terminal_stop_codons=False,
                                                        translate_target_file=False,
                                                        logger=None):
    """
    Takes a nucleotide target file and checks for unexpected stop codons when seqs are translated in the first
    forward frame. Also checks whether a seq is a multiple of three (i.e. whole codons only).

    :param str targetfile: targetfile filename
    :param bool no_terminal_stop_codons: if True, do not allow a single stop codon at C-terminal end of protein
    :param bool translate_target_file: if True, nucleotide targetfile will be translated
    :param logging.Logger logger: a logger object
    :return:
    """

    with open(targetfile, 'r') as target_file_handle:
        seqs = list(SeqIO.parse(target_file_handle, 'fasta'))

    translated_seqs_to_write = []
    seqs_needed_padding_dict = defaultdict(list)
    seqs_with_stop_codons_dict = defaultdict(list)
    seqs_with_terminal_stop_codon_dict = defaultdict(list)

    if translate_target_file:
        for seq in seqs:
            gene_name = seq.name.split('-')[-1]
            sequence, needed_padding = pad_seq(seq)
            translated_seq = sequence.seq.translate()
            if translate_target_file:
                record = SeqRecord.SeqRecord(translated_seq, id=seq.id, description='')
                translated_seqs_to_write.append(record)
            num_stop_codons = translated_seq.count('*')

            if needed_padding:
                seqs_needed_padding_dict[gene_name].append(seq)

            if num_stop_codons == 0 or \
                    (num_stop_codons == 1 and
                     re.search('[*]', str(translated_seq)[-1]) and not
                     no_terminal_stop_codons):
                seqs_with_terminal_stop_codon_dict[gene_name].append(seq)

            elif num_stop_codons >= 1:
                seqs_with_stop_codons_dict[gene_name].append(seq)

    return (seqs_with_terminal_stop_codon_dict,
            translated_seqs_to_write,
            seqs_with_stop_codons_dict,
            seqs_needed_padding_dict)


def check_targetfile(targetfile,
                     targetfile_type,
                     full_sample_directory=None,
                     using_bwa=False,
                     skip_targetfile_checks=False,
                     no_terminal_stop_codons=False,
                     sliding_window_size=0,
                     complexity_minimum_threshold=0.0,
                     running_as_subcommand=False,
                     logger=None):
    """
    - Checks target-file fasta header formatting ("taxon*-unique_gene_ID").
    - Checks for duplicate gene names in the targetfile.
    - Reports the number of unique genes (each can have multiple representatives) in the targetfile.
    - Checks that seqs in target file can be translated from the first codon position in the forwards frame (multiple of
      three, no unexpected stop codons), and logs a warning if not.
    - If not running_as_subcommand: if targetfile is DNA but using_bwa is False, translate the targetfile, write it
      to the sample directory with name 'translated_target_file.fasta', and return the path

    :param str targetfile: full path to the targetfile
    :param str targetfile_type: string describing target file sequence type i.e 'DNA' or 'protein'
    :param bool using_bwa: True if the --bwa flag is used; a nucleotide target file is expected in this case
    :param bool skip_targetfile_checks: if True, skip checks (but translate DNA file if necessary)
    :param bool no_terminal_stop_codons: when testing for open reading frames, do not allow a translated frame to have
     a single stop codon at the C-terminus of the translated protein sequence. Default is False. Applies to '
     hybpiper check_targetfile' only
    :param int sliding_window_size: number of characters (single-letter DNA or amino-acid codes) to include in the
     sliding window when checking for sequences with  low-complexity-regions. Applies to 'hybpiper check_targetfile'
     only
    :param float complexity_minimum_threshold: minimum threshold value. Beneath this value, the sequence in the
     sliding window is flagged as low complexity, and the corresponding target file sequence is reported as having
     low-complexity regions. Applies to 'hybpiper check_targetfile' only
    :param path full_sample_directory: path to the sample directory
    :param bool running_as_subcommand: True if running 'hybpiper check_targetfile' only
    :param logging.Logger logger: a logger object
    :return: None, str: NoneType or path to the targetfile (translated if necessary)
    """

    # Set defaults:
    error_with_targetfile = False
    potential_issue_with_targetfile = False
    translate_target_file_write = False
    low_complexity_sequences = None

    # Determine whether a translation check is needed:
    if targetfile_type == 'DNA':
        translate_target_file = True
    elif targetfile_type == 'protein':
        translate_target_file = False

    if running_as_subcommand:
        print(f'{"[INFO]:":10} HybPiper version {__version__} was called with these arguments:')
        fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                             break_on_hyphens=False)
        print(f'{fill}\n')

    else:
        if using_bwa and targetfile_type == 'protein':
            fill = textwrap.fill(f'{"[ERROR]:":10} You have specified that your target file contains protein '
                                 f'sequences but provided the flag --bwa. You need a nucleotide target file to use BWA '
                                 f'for read mapping!',
                                 width=90, subsequent_indent=" " * 11)
            log_or_print(fill, logger=logger, logger_level='error')
            sys.exit(1)

        elif not using_bwa and targetfile_type == 'DNA':
            fill = textwrap.fill(f'{"[WARNING]:":10} You have specified that your target file contains DNA '
                                 f'sequences, but BLASTx or DIAMOND has been selected for read mapping. Your '
                                 f'target file will be translated!', width=90, subsequent_indent=' ' * 11)
            logger.info(f'{fill}')
            translate_target_file_write = True

    if skip_targetfile_checks:
        log_or_print(f'{"[WARNING]:":10} Skipping target file checks!', logger=logger,
                     logger_level='warning')

        if translate_target_file_write:  # Can only be True when not running as subcommand
            translated_target_file = f'{full_sample_directory}/translated_target_file.fasta'
            fill = fill_forward_slash(f'{"[INFO]:":10} Writing a translated target file to:'
                                      f' {translated_target_file}',
                                      width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                      break_on_forward_slash=True)
            logger.info(f'{fill}')

            translated_seqs_to_write = []
            with open(targetfile, 'r') as target_file_handle:
                seqs = list(SeqIO.parse(target_file_handle, 'fasta'))

                for seq in seqs:
                    sequence, needed_padding = pad_seq(seq)
                    translated_seq = sequence.seq.translate()
                    record = SeqRecord.SeqRecord(translated_seq, id=seq.id, description='')
                    translated_seqs_to_write.append(record)

            # Write the translated target file to the sample directory:
            with open(f'{translated_target_file}', 'w') as translated_handle:
                SeqIO.write(translated_seqs_to_write, translated_handle, 'fasta')

            return translated_target_file  # i.e. use translated file path for return value
        else:
            return targetfile
    else:
        log_or_print(f'{"[INFO]:":10} Checking target file for issues...', logger=logger)

    # Set directory and file name for targetfile report:
    parent_dir = os.path.dirname(targetfile)
    target_file_basename = os.path.basename(targetfile)
    fn, ext = os.path.splitext(target_file_basename)

    if running_as_subcommand:  # i.e. no sample directory - write report to current working directory
        cwd = os.getcwd()
        outfile_check_report = f'{cwd}/check_targetfile_report-{fn}.txt'
    else:
        outfile_check_report = f'{full_sample_directory}/check_targetfile_report-{fn}.txt'

    # Check that seqs in target file can be translated from the first codon position in the forwards frame:
    (seqs_with_terminal_stop_codon_dict,
     translated_seqs_to_write,
     seqs_with_stop_codons_dict,
     seqs_needed_padding_dict) = \
        check_target_file_stop_codons_and_multiple_of_three(targetfile,
                                                            no_terminal_stop_codons=no_terminal_stop_codons,
                                                            translate_target_file=translate_target_file,
                                                            logger=logger)

    if translate_target_file_write:  # Can only be True when not running as subcommand
        translated_target_file = f'{full_sample_directory}/translated_target_file.fasta'
        fill = fill_forward_slash(f'{"[INFO]:":10} Writing a translated target file to:'
                                  f' {translated_target_file}',
                                  width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                  break_on_forward_slash=True)
        logger.info(f'{fill}')

        with open(f'{translated_target_file}', 'w') as translated_handle:
            SeqIO.write(translated_seqs_to_write, translated_handle, 'fasta')

        targetfile_to_return = translated_target_file # i.e. use translated file path for return value
    else:
        targetfile_to_return = targetfile  # input target file path

    # Check targetfile header and duplicate gene names:
    (incorrectly_formatted_fasta_headers,
     duplicated_gene_names,
     gene_lists) = check_target_file_headers_and_duplicate_names(targetfile)

    # Check  for low-complexity sequences only when running check_targetfile as a subcommand:
    if running_as_subcommand:
        if targetfile_type == 'DNA':
            fill = textwrap.fill(f'{"[INFO]:":10} The target file {targetfile} has been specified as containing '
                                 f'DNA sequences. These DNA sequences will be checked for low-complexity regions. '
                                 f'NOTE: the sequences flagged as having low-complexity regions can sometimes differ '
                                 f'between a DNA target file and the corresponding translated protein target file. If '
                                 f'you translate your target file to run "hybpiper assemble" with BLASTx/DIAMOND, we '
                                 f'recommend re-checking the translated sequences for low-complexity regions.',
                                 width=90,
                                 subsequent_indent=" " * 11)

            log_or_print(fill, logger=logger)

            (low_complexity_sequences,
             window_size,
             entropy_value) = \
                low_complexity_check(targetfile,
                                     targetfile_type,
                                     translate_target_file=False,
                                     window_size=sliding_window_size,
                                     entropy_value=complexity_minimum_threshold,
                                     logger=logger)

        elif targetfile_type == 'protein':
            fill = textwrap.fill(f'{"[INFO]:":10} The target file {targetfile} has been specified as containing '
                                 f'protein sequences. These protein sequences will be checked for low-complexity '
                                 f'regions',
                                 width=90, subsequent_indent=" " * 11)

            log_or_print(fill, logger=logger)

            (low_complexity_sequences,
             window_size,
             entropy_value) = \
                low_complexity_check(targetfile,
                                     targetfile_type,
                                     translate_target_file=False,
                                     window_size=sliding_window_size,
                                     entropy_value=complexity_minimum_threshold,
                                     logger=logger)

    ####################################################################################################################
    # Check for issues, log and write report:
    ####################################################################################################################
    with open(outfile_check_report, 'w') as report_handle:

        # Fasta header formatting:
        report_handle.write('[Sequences with incorrectly formatted fasta headers]\n\n')
        if len(incorrectly_formatted_fasta_headers) != 0:
            error_with_targetfile = True

            fill = textwrap.fill(f'{"[ERROR!]:":10} Some sequences in your target file have incorrectly formatted '
                                 f'fasta headers. Sequence names have been written to a report file.',
                                 width=90, subsequent_indent=' ' * 11)
            log_or_print(fill, logger=logger, logger_level='error')
            fill = textwrap.fill(f'Please see target file formatting requirements here: '
                                 f'https://github.com/mossmatters/HybPiper/wiki#12-target-file')

            log_or_print(textwrap.indent(fill, ' ' * 11), logger=logger, logger_level='error')
            log_or_print('', logger=logger, logger_level='error')

            # Write to file:
            for seq_name in incorrectly_formatted_fasta_headers:
                report_handle.write(f'{seq_name}\n')
            report_handle.write(f'\nPlease see target file formatting requirements here:\n'
                                f'https://github.com/mossmatters/HybPiper/wiki#12-target-file\n\n')
        else:
            report_handle.write(f'None\n\n')

        # Duplicated sequence names:
        report_handle.write('[Sequences with duplicated fasta headers]\n\n')
        if len(duplicated_gene_names) != 0:
            error_with_targetfile = True

            fill = textwrap.fill(f'{"[ERROR!]:":10} Some sequence names in your target file occur more than once. '
                                 f'Please remove duplicate genes before running HybPiper! Sequence names have been '
                                 f'written to a report file.',
                                 width=90, subsequent_indent=' ' * 11)

            log_or_print(fill, logger=logger, logger_level='error')

            # Write to file:
            for seq_name in duplicated_gene_names:
                report_handle.write(f'{seq_name}\n')
            report_handle.write(f'\n')
        else:
            report_handle.write(f'None\n\n')

        # Sequences with terminal stop codons:
        report_handle.write('[Sequences with terminal stop codons]\n\n')
        if seqs_with_terminal_stop_codon_dict:
            seq_list = [seq.name for gene_name, target_file_sequence_list in
                        seqs_with_terminal_stop_codon_dict.items() for seq in target_file_sequence_list]
            fill = textwrap.fill(
                f'{"[INFO]:":10} There are {len(seq_list)} sequences in your target file that contain a single '
                f'terminal stop codon. Sequence names have been written to a report file. ',
                width=90, subsequent_indent=' ' * 11)

            log_or_print(f'{fill}', logger=logger, logger_level='info')

            # Write to file:
            for seq_name in seq_list:
                report_handle.write(f'{seq_name}\n')
            report_handle.write(f'\n')
        else:
            report_handle.write(f'None\n\n')

        # Sequences with unexpected stop codons:
        report_handle.write('[Sequences with unexpected stop codons]\n\n')
        if seqs_with_stop_codons_dict:
            potential_issue_with_targetfile = True

            seq_list = [seq.name for gene_name, target_file_sequence_list in seqs_with_stop_codons_dict.items() for seq
                        in target_file_sequence_list]
            fill = textwrap.fill(
                f'{"[WARNING]:":10} There are {len(seq_list)} sequences in your target file that contain unexpected '
                f'stop codons when translated in the first forwards frame. If your target file contains only '
                f'protein-coding sequences, please check these sequences, and/or run "hybpiper fix_targetfile". '
                f'Sequence names have been written to a report file.',
                width=90, subsequent_indent=' ' * 11)

            log_or_print(f'{fill}', logger=logger)

            # Write to file:
            for seq_name in seq_list:
                report_handle.write(f'{seq_name}\n')
            report_handle.write(f'\n')
        else:
            report_handle.write(f'None\n\n')

        # Sequences not a multiple of three:
        report_handle.write('[Sequences that are not a multiple of three]\n\n')
        if seqs_needed_padding_dict:
            seq_list = [seq.name for gene_name, target_file_sequence_list in seqs_needed_padding_dict.items() for seq
                        in target_file_sequence_list]
            fill = textwrap.fill(
                f'{"[WARNING]:":10} There are {len(seq_list)} sequences in your target file that are not multiples of '
                f'three. If your target file contains only protein-coding sequences, please check these sequences, '
                f'and/or run "hybpiper fix_targetfile". Sequence names have been written to a report file.',
                width=90, subsequent_indent=' ' * 11)

            log_or_print(f'{fill}', logger=logger)

            # Write to file:
            for seq_name in seq_list:
                report_handle.write(f'{seq_name}\n')
            report_handle.write(f'\n')
        else:
            report_handle.write(f'None\n\n')

        # Low complexity sequences:
        if running_as_subcommand:
            report_handle.write('[Sequences with low complexity regions]\n\n')

            if low_complexity_sequences:
                potential_issue_with_targetfile = True

                fill_1 = textwrap.fill(
                    f'{"[WARNING]:":10} The target file provided ({os.path.basename(targetfile)}) contains '
                    f'sequences with low-complexity regions. Sequence names have been written to a report file. These '
                    f'sequences can cause problems when running HybPiper, see:',
                    width=90, subsequent_indent=" " * 11)

                fill_2 = (f'{" " * 10} https://github.com/mossmatters/HybPiper/wiki/Troubleshooting,-common-issues,'
                          f'-and-recommendations#12-target-files-with-low-complexity-sequences-troubleshooting')

                fill_3 = f'{" " * 10} We recommend one of the following approaches:'

                fill_4 = textwrap.fill(
                    f'1) Remove these sequence from your target file, ensuring that your file still contains other '
                    f'representative sequences for the corresponding genes. This can be done manually, or via the '
                    f'command "hybpiper fix_targetfile", see',
                    width=90, initial_indent=" " * 11, subsequent_indent=" " * 14)

                fill_5 = (f'{" " * 13} https://github.com/mossmatters/HybPiper/wiki/Troubleshooting,-common-issues,'
                          f'-and-recommendations#14-fixing-and-filtering-your-target-file')

                fill_6 = textwrap.fill(
                    f'2) Start the run using the parameter "--timeout_assemble" (e.g. "--timeout_assemble '
                    f'200"). See details at:',
                    width=90, initial_indent=" " * 11, subsequent_indent=" " * 14, break_on_hyphens=False)

                fill_7 = (f'{" " * 13} https://github.com/mossmatters/HybPiper/wiki/Full-pipeline-parameters#10'
                          f'-hybpiper-assemble')

                log_or_print(f'{fill_1}\n{fill_2}\n\n{fill_3}\n\n{fill_4}\n{fill_5}\n\n{fill_6}\n{fill_7}\n')

                # Write to file:
                for seq_name in low_complexity_sequences:
                    report_handle.write(f'{seq_name}\n')
                report_handle.write(f'\n')
            else:
                report_handle.write(f'None\n\n')

        else:
            report_handle.write('[Sequences with low complexity regions]\n\n')
            report_handle.write('Check not performed. Please run `hybpiper check_targetfile` if you would like to '
                                'perform this check (recommended).\n')

    if not error_with_targetfile and not potential_issue_with_targetfile:
        log_or_print(f'{"[INFO]:":10} No problems found in target file!', logger=logger)

    log_or_print(f'{"[INFO]:":10} A target file report has been written to:', logger=logger)
    log_or_print(f'{" " * 10} {outfile_check_report}')

    if running_as_subcommand:
        # Write a control file with current settings and any low-complexity sequence names; used as input to `hybpiper
        # fix_targetfile`:
        write_fix_targetfile_controlfile(targetfile_type,
                                         translate_target_file,
                                         no_terminal_stop_codons,
                                         low_complexity_sequences,
                                         window_size,
                                         entropy_value)

    if error_with_targetfile:
        sys.exit(1)

    return targetfile_to_return


def check_for_previous_run_output(full_sample_directory,
                                  start_from,
                                  end_with,
                                  assemble_stages_dict):
    """
    Checks for output files from a previous run.

    :param path full_sample_directory: full path to the sample directory
    :param str start_from: pipeline step to start from
    :param str end_with: pipeline step to end with
    :param dict assemble_stages_dict: dict of assembly stages and corresponding integer

    :return dict previous_files_dict: dictionary of previous files for each pipeline step >= start_from and <= end_with
    """

    previous_files_dict = defaultdict(list)

    # Search for files from step map_reads:
    bam = glob.glob(f'{full_sample_directory}/*.bam')
    blastx = glob.glob(f'{full_sample_directory}/*.blastx')

    # Search for files from step distribute_reads:
    pre_existing_fastas = glob.glob(f'{full_sample_directory}/*/*_interleaved.fasta') + \
                          glob.glob(f'{full_sample_directory}/*/*_unpaired.fasta')

    # Search for files from step assemble_reads:
    pre_existing_assemblies = glob.glob(f'{full_sample_directory}/*/*_contigs.fasta')

    # Search for files from step extract_contigs:
    fna_sequences = glob.glob(f'{full_sample_directory}/*/*/sequences/FNA/*.FNA')

    # Populate dictionary with any files found for each pipeline step:
    if bam:
        previous_files_dict['map_reads'].append(bam)
    if blastx:
        previous_files_dict['map_reads'].append(blastx)
    if pre_existing_fastas:
        previous_files_dict['distribute_reads'].append(pre_existing_fastas)
    if pre_existing_assemblies:
        previous_files_dict['assemble_reads'].append(pre_existing_assemblies)
    if fna_sequences:
        previous_files_dict['extract_contigs'].append(pre_existing_assemblies)

    # Filter the dict to retain only selected pipeline steps:
    start_from_int = assemble_stages_dict[start_from]
    end_with_int = assemble_stages_dict[end_with]

    previous_files_dict_filtered = copy.deepcopy(previous_files_dict)

    if start_from_int == end_with_int:  # i.e. only running a single step
        for key in previous_files_dict.keys():
            if assemble_stages_dict[key] != start_from_int:
                del previous_files_dict_filtered[key]
    else:
        for key in previous_files_dict.keys():
            if assemble_stages_dict[key] not in range(start_from_int, end_with_int +1):
                del previous_files_dict_filtered[key]

    return previous_files_dict_filtered


def parse_compressed_sample(sample_name,
                            sampledir_parent,
                            lock,
                            counter):
    """
    Parses a *.tar.gz compressed tarball sample directory and return a list of all folder and files within the tarball.

    :param str sample_name: name of the sample
    :param str sampledir_parent: name of the sample parent directory
    :param multiprocessing.managers.AcquirerProxy lock:
    :param multiprocessing.managers.ValueProxy counter:
    :return list member_names: a list of all contents of the <sample_id>.tar.gz
    """

    compressed_sample = f'{sampledir_parent}/{sample_name}.tar.gz'

    with tarfile.open(compressed_sample, 'r:gz') as tarfile_handle:

        member_name_and_size_dict = dict()

        for tarinfo in tarfile_handle.getmembers():
            member_name_and_size_dict[tarinfo.name] = tarinfo.size

    with lock:
        counter.value += 1
        return (sample_name,
                member_name_and_size_dict,
                counter.value)


def get_compressed_seqrecords(tarfile_handle,
                              member_name):
    """
    Returns a list of seqrecords from fasta sequences within a multi-fasta file located within a compressed tarball of
    a given sample directory.

    :param tarfile_handle:
    :param member_name: name of the tarfile member to find
    :return list: returns a list of seqrecords
    """

    extracted_file = tarfile_handle.extractfile(member_name)
    lines = extracted_file.read().decode('utf-8', errors='ignore')
    seq_io = io.StringIO(lines)
    seqrecords = list(SeqIO.parse(seq_io, 'fasta'))

    return seqrecords


def get_compressed_file_lines(tarfile_handle,
                              member_name):
    """

    :param tarfile_handle:
    :param member_name: name of the tarfile member to extract
    :return:
    """

    extracted_file = tarfile_handle.extractfile(member_name)
    lines = extracted_file.read().decode('utf-8', errors='ignore')
    lines = lines.split('\n')
    lines = [line for line in lines if line]

    return lines


def get_compressed_file_lines_generator(tarfile_handle,
                                        member_name):
    """

    :param tarfile_handle:
    :param member_name: name of the tarfile member to extract
    :return:
    """
    with tarfile_handle.extractfile(member_name) as f:
        # Read and yield lines one by one
        for line in f:
            yield line.decode('utf-8', errors='ignore').rstrip()


def get_bamtools_flagstat_lines_from_compressed(sample_name,
                                                sampledir_parent,
                                                bam_file):
    """

    :param sample_name: sample name
    :param sampledir_parent: sampledir_parent name
    :param bam_file: name of the bam file to extract to disk
    :return:
    """

    # Create temp output directory:
    outdir = f'{os.getcwd()}/temp_bam_files'

    message_list = [f'No samtools flagstat output for file "{bam_file}" was found for sample {sample_name}. '
                    f'Attempting to extract {sample_name}.bam file to directory "{outdir}"...']

    try:
        outdir = createfolder(outdir)
    except OSError:
        message_list.append(f'Could not create directory {outdir}. Do you have write permission for the parent '
                            f'directory?')
        return ([],
                message_list)

    # Extract the bam file to a temp directory:
    with tarfile.open(f'{sampledir_parent}/{sample_name}.tar.gz', 'r:gz') as tarfile_handle:
        tarfile_handle.extract(bam_file, outdir)

    message_list.append(f'Success! Running `samtools flagstat` now...')

    # Run samtools flagstat:
    samtools_cmd = f'samtools flagstat --output-fmt tsv {outdir}/{bam_file}'

    try:
        result = subprocess.run(samtools_cmd,
                                shell=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                universal_newlines=True,
                                check=True)

        message_list.append(f'Success!')

    except subprocess.CalledProcessError as exc:
        message_list.append(f'samtools flagstat mapping FAILED. Output is: {exc}')
        message_list.append(f'samtools flagstat mapping stdout is: {exc.stdout}')
        message_list.append(f'samtools flagstat mapping stderr is: {exc.stderr}')
        return ([],
                message_list)

    # Remove the temp directory:
    shutil.rmtree(outdir)

    return ([line.rstrip() for line in result.stdout.split('\n')],
            message_list)


def get_bamtools_flagstat_lines_from_uncompressed(sample_name,
                                                  sampledir_parent,
                                                  bam_file):
    """

    :param sample_name: sample name
    :param sampledir_parent: sampledir_parent name
    :param bam_file: name of the bamfile to run `samtools flagstat` on
    :return:
    """

    message_list = [f'No samtools flagstat output for file "{bam_file}" was found for sample {sample_name}. Running '
                    f'`samtools flagstat` now...']

    # Run samtools flagstat:
    samtools_cmd = f'samtools flagstat --output-fmt tsv {sampledir_parent}/{bam_file}'

    try:
        result = subprocess.run(samtools_cmd,
                                shell=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                universal_newlines=True,
                                check=True)

        message_list.append(f'Success!')

    except subprocess.CalledProcessError as exc:
        message_list.append(f'samtools flagstat mapping FAILED. Output is: {exc}')
        message_list.append(f'samtools flagstat mapping stdout is: {exc.stdout}')
        message_list.append(f'samtools flagstat mapping stderr is: {exc.stderr}')
        return ([],
                message_list)

    return ([line.rstrip() for line in result.stdout.split('\n')],
            message_list)


def check_namelist(namelist,
                   logger):
    """
    Parse namelist text file and check for issues

    :param namelist:
    :param logger:
    :return:
    """

    set_of_sample_names = set()
    duplicate_sample_names_set = set()
    samples_with_slashes_list = []

    with open(namelist, 'r') as namelist_handle:
        for line in namelist_handle.readlines():
            sample_name = line.rstrip()
            if sample_name:
                if sample_name in set_of_sample_names:
                    duplicate_sample_names_set.add(sample_name)
                if re.search('/', sample_name):
                    samples_with_slashes_list.append(sample_name)

                set_of_sample_names.add(sample_name)

    if len(samples_with_slashes_list) != 0:
        logger.error(f'{"[ERROR]:":10} A sample name must not contain forward slashes. The file "{namelist}" '
                     f'contains the following sample names with forward slashes:\n')
        for sample_name in sorted(samples_with_slashes_list):
            logger.error(f'{" " * 10} {sample_name}')
        logger.error(f'')

    if len(duplicate_sample_names_set) != 0:
        logger.error(f'{"[ERROR]:":10} The file "{namelist}" should not contain duplicate sample names. '
                     f'More than one entry was found for the following sample names:\n')
        for sample_name in sorted(duplicate_sample_names_set):
            logger.error(f'{" " * 10} {sample_name}')
        logger.error(f'')

    if len(samples_with_slashes_list) != 0 or len(duplicate_sample_names_set) != 0:
        sys.exit()

    return set_of_sample_names


def check_for_compressed_and_uncompressed_samples(set_of_sample_names,
                                                  sampledir_parent,
                                                  logger):
    """
    Check if there is both an uncompressed folder AND a compressed file for any sample

    :param set_of_sample_names:
    :param sampledir_parent: parent directory containing sample files/folders
    :param logger:
    :return:
    """

    ####################################################################################################################
    samples_found = []
    compressed_samples_set = set()
    uncompressed_samples_set = set()

    for sample_name in set_of_sample_names:

        compressed_sample = f'{sampledir_parent}/{sample_name}.tar.gz'
        uncompressed_sample = f'{sampledir_parent}/{sample_name}'

        if os.path.isfile(compressed_sample):
            compressed_samples_set.add(sample_name)
            samples_found.append(sample_name)

        if os.path.isdir(uncompressed_sample):
            uncompressed_samples_set.add(sample_name)
            samples_found.append(sample_name)

    both_compressed_and_uncompressed_present = compressed_samples_set.intersection(uncompressed_samples_set)

    if len(both_compressed_and_uncompressed_present) != 0:
        logger.error(f'{"[ERROR]:":10} Both a compressed and an un-compressed sample folder have been found for the '
                     f'following samples - please remove one:\n')
        for sample_name in sorted(list(both_compressed_and_uncompressed_present)):
            logger.error(f'{" " * 10} {sample_name}')
        logger.error(f'')
        sys.exit()

    return (samples_found,
            compressed_samples_set,
            uncompressed_samples_set)


def check_for_missing_samples(namelist,
                              set_of_sample_names,
                              samples_found,
                              sampledir_parent,
                              logger):
    """
    Check for sample names present in the namelist.txt but not found in the directory provided

    :param namelist:
    :param set_of_sample_names:
    :param samples_found:
    :param sampledir_parent: parent directory containing sample files/folders
    :param logger:
    :return:
    """

    samples_missing = sorted(list(set(set_of_sample_names) - set(samples_found)))

    if samples_missing:

        fill = fill_forward_slash(f'{"[WARNING]:":10} File {namelist} contains samples not found in directory '
                                  f'"{sampledir_parent}". The missing samples are:',
                                  width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                  break_on_forward_slash=True)

        logger.warning(f'{fill}\n')

        for name in samples_missing:
            logger.warning(f'{" " * 10} {name}')
        logger.warning('')

    list_of_sample_names = sorted([x for x in set_of_sample_names if x not in samples_missing])

    return list_of_sample_names
