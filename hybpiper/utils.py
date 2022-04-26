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
            logger.debug(string)
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
    :return set low_entropy_seqs: a set of sequences that contain low entropy substrings
    """

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

    low_entropy_seqs = set()
    for seq in SeqIO.parse(targetfile, "fasta"):
        for i in range(0, len(seq.seq) - (window_size - 1)):
            window_seq = str(seq.seq[i:i + window_size])
            window_shannon_entropy = shannon_entropy(window_seq)
            # print(window_seq, window_shannon_entropy)
            if window_shannon_entropy <= entropy_value:
                low_entropy_seqs.add(seq.name)
                break

    return low_entropy_seqs
