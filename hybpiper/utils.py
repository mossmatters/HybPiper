#!/usr/bin/env python

"""
This module contains some general functions and classes.
"""

import re
from textwrap import TextWrapper
import os
import collections
import scipy
import textwrap
from Bio import SeqIO
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
    # for seq in SeqIO.parse(targetfile, "fasta"):
    #     total_number_of_seqs += 1

    low_entropy_seqs = set()
    for sequence in progressbar.progressbar(total_seqs, max_value=len(total_seqs), min_poll_interval=30,
                                            widgets=widgets):

    # for seq in SeqIO.parse(targetfile, "fasta"):
        for i in range(0, len(sequence.seq) - (window_size - 1)):
            window_seq = str(sequence.seq[i:i + window_size])
            window_shannon_entropy = shannon_entropy(window_seq)
            if window_shannon_entropy <= entropy_value:
                low_entropy_seqs.add(sequence.name)
                break

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

                    log_or_print(f'\n{exe} version 2.4.0 is required, but your version is {version}. Please update '
                                 f'your Exonerate version!\n', logger=logger)
        else:
            log_or_print(f'{exe:20} not found in your $PATH!', logger=logger)
            everything_is_awesome = False

    log_or_print('', logger=logger)

    return everything_is_awesome


def make_basename(readfiles, prefix=None):
    """
    Unless prefix is set, generate a directory based on the readfiles name, using everything up to the first
    underscore. If prefix is set, generate the directory "prefix" and set basename to be the last component of the path.

    :param list readfiles: one or more read files used as input to the pipeline
    :param str prefix: directory name for sample pipeline output
    :return str parent directory, directory name
    """

    if prefix:
        if not os.path.exists(prefix):
            os.makedirs(prefix)
        prefixparendir, prefix = os.path.split(prefix)
        if not prefix:
            # if prefix has a trailing /, prefixparendir will have the / stripped and prefix will be empty,
            # so try again
            prefix = os.path.split(prefixparendir)[1]
        return prefixparendir, prefix

    # --prefix is not set on cmd line;  Write output to subdir in "."
    basename = os.path.split(readfiles[0])[1].split('_')[0]
    if not os.path.exists(basename):
        os.makedirs(basename)
    return '.', basename


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


def createfolder(directory):
    """
    Attempts to create a directory named after the name provided, and provides an error message on failure
    """
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print(f'Error: Creating directory: {directory}')


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

