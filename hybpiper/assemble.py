#!/usr/bin/env python

"""
 _    _            _       _____
| |  | |          | |     |  _  \
| |__| | __    __ | |___  | |_| |  _   _____   _____   _____
|  __  | \ \  / / |  __ \ |  ___/ | | |  _  \ |  .  | |  _  \
| |  | |  \ \/ /  | |__ | | |     | | | |_| | |  __/  | |
|_|  |_|   \  /   |_____/ |_|     |_| |  ___/ |_____| |_|
           / /                        | |
          /_/                         |_|

HybPiper Version 2.0 release candidate (March 2022)

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

To view available parameters and help for any subcommand, simply type e.g. 'hybpiper assemble -h'.

==> NOTE <==
The command/script 'read_first.py' no longer exists, and has been replaced by the subcommand 'assemble'. So,
if you had previously run 'reads_first.py' on a sample using the command e.g.:

    python /<path_to>/reads_first.py -t test_targets.fasta -r NZ281_R*_test.fastq --prefix NZ281 --bwa

...this is now replaced by the command:

    hybpiper assemble -t test_targets.fasta -r NZ281_R*_test.fastq --prefix NZ281 --bwa

==> NOTE <==
The recovery of introns and supercontigs, previously achieved via the script 'intronerate.py',
is now incorporated in to the 'hybpiper assemble' command. It can be enabled using the flag
'--run_intronerate', e.g.:

    hybpiper assemble -t test_targets.fasta -r NZ281_R*_test.fastq --prefix NZ281 --bwa --run_intronerate

==> NOTE <==
The command/script 'get_seq_lengths.py' no longer exists, and this functionality has been incorporated in to
the command 'hybpiper stats'. The sequence length details that were previously printed to screen are now written to
the file 'seq_lengths.tsv', by default. Similarity, the stats details that were previously written to screen by
'hybpiper_stats.py' are now written to the file 'hybpiper_stats.tsv', by default.

For full details of all commands and changes, please reads the Wiki page at **LINK** and the changelog at **LINK**.

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
from concurrent.futures.process import ProcessPoolExecutor
import multiprocessing
from multiprocessing import Manager
from concurrent.futures import wait, as_completed
import pkg_resources

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
from hybpiper import distribute_reads_to_targets_bwa
from hybpiper import distribute_reads_to_targets
from hybpiper import distribute_targets
from hybpiper import spades_runner
from hybpiper import exonerate_hits
from hybpiper import hybpiper_stats
from hybpiper import retrieve_sequences
from hybpiper import paralog_retriever
from hybpiper import gene_recovery_heatmap
from hybpiper import hybpiper_subparsers


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


def seq_is_dna_or_protein(seq, alphabet='dna', ambiguity_characters=None):
    """
    Check that a sequence only contains values from the DNA alphabet.

    Adapted from https://www.biostars.org/p/102/

    :param str seq: a string of letters representing a DNA or amino-acid sequence
    :param str alphabet: the type of alphabet to search seq for; choices are 'dna' or 'protein'
    :param str ambiguity_characters: a string of allowed ambiguity characters
    :return bool: True if given seq only contains values from the DNA alphabet
    """

    if alphabet not in ['dna', 'protein']:
        raise ValueError(f'alphabet should be either "dna" or "protein". Current value is: {alphabet} ')

    dna_pattern = f'^[acgtn{ambiguity_characters}]*$'
    protein_pattern = f'^[acdefghiklmnpqrstvwyx{ambiguity_characters}]*$'

    alphabets = {'dna': re.compile(dna_pattern, re.IGNORECASE),
                 'protein': re.compile(protein_pattern, re.IGNORECASE)}

    if alphabets[alphabet].search(seq) is not None:
        return True
    else:
        return False


def check_dependencies(logger=None):
    """
    Checks for the presence of executables and Python packages. Returns a boolean.

    :param logging.Logger logger: a logger object
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
                   'diamond']

    if not logger:
        print(f'{"[NOTE]:":10} Checking for external dependencies:\n')
    else:
        logger.info(f'{"[NOTE]:":10} Checking for external dependencies:\n')

    everything_is_awesome = True
    for exe in executables:
        exe_loc = py_which(exe)
        if exe_loc:
            if not logger:
                print(f'{exe:20} found at {exe_loc}')
            else:
                logger.info(f'{exe:20} found at {exe_loc}')
        else:
            if not logger:
                print(f'{exe:20} found at {exe_loc}')
            else:
                logger.info(f'{exe:20} not found in your $PATH!')
            everything_is_awesome = False
    if not logger:
        print('')
    else:
        logger.info('')
    return everything_is_awesome


def check_targetfile(targetfile, using_bwa, targetfile_ambiguity_codes, logger=None):
    """
    - Checks target-file fasta header formatting ("taxon*-unique_gene_ID").
    - Reports the number of unique genes (each can have multiple representatives) in the targetfile.
    - Performs a detection of whether targetfile is DNA or amino-acid.
    - Checks that seqs in target file can be translated from the first codon position in the forwards frame (multiple of
      three, no unexpected stop codons), and logs a warning if not.
    - If targetfile comprises nucleotides but using_bwa is False, translate the targetfile and return the path

    :param str targetfile: path to the targetfile
    :param bool using_bwa: True if the --bwa flag is used; a nucleotide target file is expected in this case
    :param str targetfile_ambiguity_codes: a string of allowed ambiguity codes when checking if targetfile is DNA
    :param logging.Logger logger: a logger object
    :return: None, str: NoneType or path to the translated targetfile
    """

    # Check target file fasta header formatting:
    logger.info(f'{"[NOTE]:":10} Checking target file FASTA header formatting...')

    target_file_path, target_file_name = os.path.split(targetfile)
    file_name, ext = os.path.splitext(target_file_name)
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
        logger.error(f'{"[ERROR!]:":10} The following sequences in your target file have incorrectly formatted fasta '
                     f'headers:\n')
        fill = textwrap.fill(f'{seq_list}')
        logger.info(textwrap.indent(fill, ' ' * 11))
        logger.info('')
        sys.exit(1)  # target file fasta header formatting should be fixed!
    else:
        logger.info(f'{"[NOTE]:":10} The target file FASTA header formatting looks good!')

    # Check for duplicated genes:
    duplicated_genes = []
    for gene, gene_count in check_for_duplicate_genes_dict.items():
        if gene_count > 1:
            duplicated_genes.append(gene)
    if duplicated_genes:
        gene_list = ' '.join(duplicated_genes)
        logger.error(f'{"[ERROR!]:":10} The following sequences in your target file occur more than once:\n')
        fill = textwrap.fill(f'{gene_list}')
        logger.info(textwrap.indent(fill, ' ' * 11))
        logger.error(f'\nPlease remove duplicate genes before running HybPiper!')
        sys.exit(1)  # duplicate genes in targetfile should be removed!

    # Report the number of unique genes represented in the target file:
    logger.info(f'{"[NOTE]:":10} The target file contains at least one sequence for {len(gene_lists)} '
                f'unique genes.')

    # Detect whether the target file is DNA or amino-acid:
    translate_target_file = False
    number_of_target_sequences = len(seqs)
    number_of_dna_sequences = len(
        [seq_is_dna_or_protein(str(seq.seq), alphabet='dna', ambiguity_characters=targetfile_ambiguity_codes)
         for seq in seqs if
         seq_is_dna_or_protein(str(seq.seq), alphabet='dna', ambiguity_characters=targetfile_ambiguity_codes)])

    logger.debug(f'There are {number_of_target_sequences} sequences in the target file, of which'
                 f' {number_of_dna_sequences} appear to be DNA')

    if number_of_dna_sequences == number_of_target_sequences:
        logger.info(f'{"[NOTE]:":10} The target file appears to contain nucleotide sequences.')
        target_file_id = 'dna'
    else:  # Check if all seqs contain only the protein alphabet (and optional ambiguity codes)
        undetermined_seq_names = []
        for seq in seqs:
            try:
                assert seq_is_dna_or_protein(
                    str(seq.seq), alphabet='protein', ambiguity_characters=targetfile_ambiguity_codes)
            except AssertionError:
                logger.debug(f'Target file sequence {seq.id} can not be identified as nucleotide or protein, check!')
                undetermined_seq_names.append(seq.id)
        if undetermined_seq_names:
            joined_names = ' '.join(undetermined_seq_names)
            sys.exit(f'{"[ERROR!]:":10} HybPiper can not determine whether the target file provided contains '
                     f'nucleotide or protein sequences. Please check your target file. If your sequences contain '
                     f'ambiguity codes, please supply them as a string using the parameter '
                     f'"--targetfile_ambiguity_codes". The undetermined sequences are:\n{joined_names}')
        else:
            logger.info(f'{"[NOTE]:":10} The target file appears to contain protein sequences.')
            target_file_id = 'protein'

    if using_bwa and target_file_id == 'protein':
        sys.exit(f'{"[ERROR]:":10} Your target file appears to contain protein sequences. You need a nucleotide '
                 f'target file for BWA!')
    elif not using_bwa and target_file_id == 'dna':
        logger.info(f'{"[WARN!]:":10} The target file appears to contain nucleotide sequences, but BLASTx or DIAMOND '
                    f'has been selected for read mapping. Translating the target file...')
        translate_target_file = True

    # Check that seqs in target file can be translated from the first codon position in the forwards frame:
    if using_bwa or translate_target_file:
        translated_seqs_to_write = []

        seqs_needed_padding_dict = defaultdict(list)
        seqs_with_stop_codons_dict = defaultdict(list)

        for seq in seqs:
            gene_name = seq.name.split('-')[-1]
            sequence, needed_padding = distribute_targets.pad_seq(seq)
            translated_seq = sequence.seq.translate()
            if translate_target_file:
                record = SeqRecord.SeqRecord(translated_seq, id=seq.id, description='')
                translated_seqs_to_write.append(record)
            num_stop_codons = translated_seq.count('*')

            if needed_padding:
                seqs_needed_padding_dict[gene_name].append(seq)

            if num_stop_codons == 1 and re.search('[*]', str(translated_seq)[-5:]):
                logger.debug(f'Translated sequence {seq.name} contains a single stop codon in the last 5 amino-acids, '
                             f'proceeding...')
            elif num_stop_codons >= 1:
                seqs_with_stop_codons_dict[gene_name].append(seq)

        if seqs_with_stop_codons_dict:
            seq_list = [seq.name for gene_name, target_file_Sequence_list in seqs_with_stop_codons_dict.items() for seq
                        in target_file_Sequence_list]
            logger.info(f'{"[WARN!]:":10} There are {len(seq_list)} sequences in your target file that contain '
                        f'unexpected stop codons when translated in the first forwards frame. \n{" " * 11}If your '
                        f'target file contains only protein-coding sequences, please check these sequences. '
                        f'\n{" " * 11}Sequence names can be found in the sample log file.\n')
            logger.debug(f'Target file sequences with unexpected stop codons: {seq_list}')

        if seqs_needed_padding_dict:
            seq_list = [seq.name for gene_name, target_file_Sequence_list in seqs_needed_padding_dict.items() for seq
                        in target_file_Sequence_list]
            logger.info(f'{"[WARN!]:":10} There are {len(seq_list)} sequences in your target file that are not '
                        f'multiples of three. \n{" " * 11}If your target file contains only protein-coding sequences, '
                        f'please check these sequences. \n{" " * 11}Sequence names can be found in the sample log '
                        f'file.\n')
            logger.debug(f'Target file sequences that are not multiples of three: {seq_list}')

        if translate_target_file:
            translated_target_file = f'{target_file_path}/{file_name}_translated{ext}'
            logger.info(f'{"[NOTE]:":10} Writing a translated target file to: {translated_target_file}')
            with open(f'{translated_target_file}', 'w') as translated_handle:
                SeqIO.write(translated_seqs_to_write, translated_handle, 'fasta')
            return translated_target_file

    return targetfile


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


def bwa(readfiles, targetfile, basename, cpu, unpaired=False, logger=None):
    """
    Conduct a BWA search of input reads against the targetfile.

    :param list readfiles: one or more read files used as input to the pipeline
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
            logger.info(f'{"[NOTE]:":10} Making nucleotide bwa index in current directory.')
            targetfiledir = os.path.split(targetfile)[0]
            if targetfiledir:
                if os.path.realpath(targetfiledir) != os.path.realpath('.'):
                    shutil.copy(targetfile, '.')
            db_file = os.path.split(targetfile)[1]
            make_bwa_index_cmd = f'bwa index {db_file}'
            logger.info(f'{"[CMD]:":10} {make_bwa_index_cmd}')

            try:
                result = subprocess.run(make_bwa_index_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True)
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

    if not cpu:
        import multiprocessing
        cpu = multiprocessing.cpu_count()  # i.e. use all cpus. Reduce by 1?

    if isinstance(readfiles, list) and len(readfiles) < 3:
        bwa_fastq = ' '.join(readfiles)
    elif isinstance(readfiles, str):
        bwa_fastq = readfiles  # i.e. a path (str) to a file of unpaired (differ\nnet to single-end) reads
    else:
        raise ValueError(f'Can not determine whether {readfiles} is single-end, paired-end or unpaired!')

    bwa_commands = ['time bwa mem', '-t', str(cpu), db_file, bwa_fastq, ' | samtools view -h -b -S - > ']
    if unpaired:
        bwa_commands.append(f'{basename}_unpaired.bam')
    else:
        bwa_commands.append(f'{basename}.bam')
    full_command = ' '.join(bwa_commands)
    logger.info(f'{"[CMD]:":10} {full_command}')

    try:
        result = subprocess.run(full_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True)
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
            logger.info(f'{"[NOTE]:":10} Making protein blastdb in current directory.')
            if os.path.split(targetfile)[0]:
                shutil.copy(targetfile, '.')
            db_file = os.path.split(targetfile)[1]
            if diamond:
                logger.info(f'{"[NOTE]:":10} Using DIAMOND instead of BLASTx!')
                if diamond_sensitivity:
                    logger.info(f'{"[NOTE]:":10} Using DIAMOND sensitivity "{diamond_sensitivity}"')
                makeblastdb_cmd = f'diamond makedb --in {db_file} --db {db_file}'
            else:
                makeblastdb_cmd = f'makeblastdb -dbtype prot -in {db_file}'
            logger.info(f'{"[NOTE]:":10} {makeblastdb_cmd}')
            try:
                result = subprocess.run(makeblastdb_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True)
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

    # Remove previous blast results if they exist (because we will be appending)
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
        if cpu:
            full_command = f"time {pipe_cmd} | parallel -j {cpu} -k --block 200K --recstart '>' --pipe " \
                           f"'{blastx_command}' >> {basename}_unpaired.blastx"
        else:
            full_command = f"time {pipe_cmd} | parallel -k --block 200K --recstart '>' --pipe '{blastx_command}' >>" \
                           f" {basename}_unpaired.blastx"
        logger.info(f'{"[CMD]:":10} {full_command}')

        try:
            result = subprocess.run(full_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True)
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
            if cpu:
                full_command = f"time {pipe_cmd} | parallel -j {cpu} -k --block 200K --recstart '>' --pipe " \
                               f"'{blastx_command}' >> {basename}.blastx"
            else:
                full_command = f"time {pipe_cmd} | parallel -k --block 200K --recstart '>' --pipe " \
                               f"'{blastx_command}' >> {basename}.blastx"

            logger.info(f'{"[CMD]:":10} {full_command}')

            try:
                result = subprocess.run(full_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True)
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
                      merged=False, logger=None):
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
    :param logging.Logger logger: a logger object
    :return: None
    """

    # Distribute reads to gene directories:
    read_hit_dict = distribute_reads_to_targets.read_sorting(blastx_outputfile)

    if len(readfiles) == 2:
        logger.info(f'{"[NOTE]:":10} In total, {len(read_hit_dict) * 2} reads from the paired-end read files '
                    f'will be distributed to gene directories')
        single_end = False
    elif len(readfiles) == 1:
        logger.info(f'{"[NOTE]:":10} In total, {len(read_hit_dict)} reads from the single-end read file will '
                    f'be distributed to gene directories')
        single_end = True
    else:
        raise ValueError(f'Can not determine whether single-end or pair-end reads were provided!')

    distribute_reads_to_targets_bwa.distribute_reads(readfiles, read_hit_dict, merged=merged,
                                                     single_end=single_end)

    if unpaired_readfile:
        up_blastx_outputfile = blastx_outputfile.replace('.blastx', '_unpaired.blastx')
        read_hit_dict_unpaired = distribute_reads_to_targets.read_sorting(up_blastx_outputfile)
        logger.info(f'{"[NOTE]:":10} In total, {len(read_hit_dict)} reads from the unpaired read file will be '
                    f'distributed to gene directories')
        distribute_reads_to_targets.distribute_reads([unpaired_readfile], read_hit_dict_unpaired,
                                                     unpaired_readfile=unpaired_readfile)

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
                   logger=None):
    """
    When using BWA mapping, distribute sample reads to their corresponding target file gene matches.

    Distribute the 'best' target file sequence (translated if necessary) to each gene directory.

    :param str bamfile: *.bam output file from BWA alignment of sample reads to the target file
    :param list readfiles: one or more read files used as input to the pipeline
    :param str targetfile: path to targetfile (i.e. the target file)
    :param str target: specific target(s) to use. Tab-delimited file (one <gene>\t<taxon_name> per line) or single
    taxon name
    :param str/bool unpaired_readfile: a path if an unpaired file has been provided, False if not
    :param str exclude: specify sequence not to be used as a target sequence for Exonerate
    :param bool merged: if True, write and distribute fastq files for merging with BBmerge.sh (in addition to fasta)
    :param logging.Logger logger: a logger object
    :return: None
    """

    # Distribute reads to gene directories:
    read_hit_dict = distribute_reads_to_targets_bwa.read_sorting(bamfile)

    if len(readfiles) == 2:
        logger.info(f'{"[NOTE]:":10} In total, {len(read_hit_dict) * 2} reads from the paired-end read files '
                    f'will be distributed to gene directories')
        single_end = False
    elif len(readfiles) == 1:
        logger.info(f'{"[NOTE]:":10} In total, {len(read_hit_dict)} reads from the single-end read file will '
                    f'be distributed to gene directories')
        single_end = True
    else:
        raise ValueError(f'Can not determine whether single-end or pair-end reads were provided!')

    distribute_reads_to_targets_bwa.distribute_reads(readfiles, read_hit_dict, merged=merged,
                                                     single_end=single_end)

    if unpaired_readfile:
        up_bamfile = bamfile.replace('.bam', '_unpaired.bam')
        read_hit_dict_unpaired = distribute_reads_to_targets_bwa.read_sorting(up_bamfile)
        logger.info(f'{"[NOTE]:":10} In total, {len(read_hit_dict)} reads from the unpaired read file will be '
                    f'distributed to gene directories')
        distribute_reads_to_targets_bwa.distribute_reads([unpaired_readfile], read_hit_dict_unpaired,
                                                         unpaired_readfile=unpaired_readfile)

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
    :param bool unpaired: True is an unpaired readfile has been provided for the sample
    :param bool merged: True if parameter --merged is used
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
    logger.debug(f'args.unpaired is:  {unpaired}')

    if unpaired:  # Create an empty unpaired file if it doesn't exist
        for gene in open('spades_genelist.txt'):
            gene = gene.rstrip()
            if os.path.isfile(f'{gene}/{gene}_interleaved.fasta'):
                if not os.path.isfile(f'{gene}/{gene}_unpaired.fasta'):
                    open(f'{gene}/{gene}_unpaired.fasta', 'a').close()

    spades_failed = spades_runner.spades_initial('spades_genelist.txt', cov_cutoff=cov_cutoff, cpu=cpu,
                                                 kvals=kvals, paired=paired, timeout=timeout, unpaired=unpaired,
                                                 merged=merged, single_cell_mode=single_cell_mode)
    logger.info(f'{"[NOTE]:":10} Finished running initial SPAdes assemblies for all genes with reads!')
    if len(spades_failed) > 0:
        with open('failed_spades.txt', 'w') as failed_spadefile:
            failed_spadefile.write('\n'.join(spades_failed))

        spades_duds = spades_runner.rerun_spades('failed_spades.txt', cov_cutoff=cov_cutoff, cpu=cpu)
        logger.info(f'{"[NOTE]:":10} Finished re-running SPAdes assemblies for genes with unsuccessful initial '
                    f'assemblies!')

        if len(spades_duds) == 0:
            logger.info(f'{"[NOTE]:":10} All SPAdes re-runs completed successfully!')

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


def done_callback(future_returned):
    """
    Callback function for ProcessPoolExecutor futures; gets called when a future is cancelled or 'done'.

    :param concurrent.futures._base.Future future_returned: future object returned by ProcessPoolExecutor
    :return: None if future cancelled or error, future_returned.result() as result if successful
    """

    if future_returned.cancelled():
        print(f'{future_returned}: cancelled')
        return future_returned.result()
    elif future_returned.done():
        error = future_returned.exception()  # returns None if no error raised
        if error:
            print(f'{future_returned}: error returned: {error}')
            return future_returned.result()
        else:
            result = future_returned.result()
            return result


def exonerate(gene_name,
              basename,
              thresh=55,
              paralog_warning_min_length_percentage=0.75,
              depth_multiplier=10,
              no_stitched_contig=False,
              bbmap_memory=1,
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
              keep_intermediate_files=False):
    """
    :param str gene_name: name of a gene that had at least one SPAdes contig
    :param str basename: directory name for sample
    :param int thresh: percent identity threshold for stitching together Exonerate results
    :param float paralog_warning_min_length_percentage: min % of a contig vs ref protein length for a paralog warning
    :param int depth_multiplier: assign long paralog as main if coverage depth <depth_multiplier> other paralogs
    :param bool no_stitched_contig: if True, don't create stitched contigs and just use longest Exonerate hit
    :param int bbmap_memory: GB memory (RAM ) to use for bbmap.sh
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
    :return: str gene_name, str prot_length OR None, None
    """

    logger = logging.getLogger()  # Assign root logger from inside the new Python process (ProcessPoolExecutor pool)
    if logger.hasHandlers():
        logger.handlers.clear()
    worker_configurer_func(gene_name)  # set up process-specific logging to file
    logger = logging.getLogger(gene_name)
    logger.setLevel(logging.DEBUG)

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

        logger.debug(f'spades_assembly_dict is: {spades_assembly_dict}')
        logger.debug(f'best_protein_ref_dict is: {best_protein_ref_dict}')

    except FileNotFoundError as e:
        logger.error(f"\n{'[ERROR!]:':10} Couldn't find an expected file for either the SPAdes assembly or the protein "
                     f"reference for gene {gene_name}, error is {e}")
        with lock:
            counter.value += 1
            sys.stderr.write(f'\r{"[NOTE]:":10} Finished running Exonerate for gene {gene_name}, {counter.value}'
                             f'/{genes_to_process}')
        return gene_name, None  # return gene_name to that log can be re-logged to main log file

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
            keep_intermediate_files=keep_intermediate_files)

        if intronerate and exonerate_result and exonerate_result.hits_filtered_by_pct_similarity_dict:
            logger.debug(f'exonerate_result.hits_subsumed_hits_removed_overlaps_trimmed_dict for gene {gene_name} is:'
                         f' {exonerate_result.hits_subsumed_hits_removed_overlaps_trimmed_dict}')
            logger.debug(f'Running intronerate')
            exonerate_hits.intronerate(exonerate_result,
                                       spades_assembly_dict,
                                       logger=logger,
                                       no_padding_supercontigs=no_padding_supercontigs,
                                       keep_intermediate_files=keep_intermediate_files)
    else:
        exonerate_result = False

    with lock:
        counter.value += 1
        sys.stderr.write(f'\r{"[NOTE]:":10} Finished running Exonerate for gene {gene_name}, {counter.value}'
                         f'/{genes_to_process}')

    if not exonerate_text_output or not exonerate_result or not exonerate_result.stitched_contig_seqrecord:
        return gene_name, None  # return gene_name to that exonerate_hits.py log can be re-logged to main log file

    return gene_name, len(exonerate_result.stitched_contig_seqrecord)


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


def exonerate_multiprocessing(genes,
                              basename,
                              thresh=55,
                              paralog_warning_min_length_percentage=0.75,
                              pool_threads=None,
                              depth_multiplier=10,
                              no_stitched_contig=False,
                              bbmap_memory=1,
                              bbmap_subfilter=7,
                              bbmap_threads=2,
                              chimeric_stitched_contig_edit_distance=5,
                              chimeric_stitched_contig_discordant_reads_cutoff=5,
                              logger=None,
                              intronerate=False,
                              no_padding_supercontigs=False,
                              keep_intermediate_files=False):
    """
    Runs the function exonerate() using multiprocessing.

    :param list genes: list of genes that had successful SPAdes runs
    :param str basename: directory name for sample
    :param int thresh: percent identity threshold for stitching together Exonerate results
    :param float paralog_warning_min_length_percentage: min % of a contig vs ref protein length for a paralog warning
    :param int pool_threads: number of threads/cpus to use for the ProcessPoolExecutor pool
    :param int depth_multiplier: assign long paralog as main if coverage depth <depth_multiplier> other paralogs
    :param bool no_stitched_contig: if True, don't create stitched contig and just use longest Exonerate hit
    :param int bbmap_memory: GB memory (RAM ) to use for bbmap.sh
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
    :return:
    """

    logger.info(f'{"[NOTE]:":10} Running exonerate_hits for {len(genes)} genes...')
    genes_to_process = len(genes)

    if not pool_threads:
        import multiprocessing
        pool_threads = multiprocessing.cpu_count()
    logger.debug(f'exonerate_multiprocessing pool_threads is: {pool_threads}')

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(exonerate,
                                      gene_name,
                                      basename,
                                      thresh=thresh,
                                      paralog_warning_min_length_percentage=paralog_warning_min_length_percentage,
                                      depth_multiplier=depth_multiplier,
                                      no_stitched_contig=no_stitched_contig,
                                      bbmap_memory=bbmap_memory,
                                      bbmap_subfilter=bbmap_subfilter,
                                      bbmap_threads=bbmap_threads,
                                      chimeric_stitched_contig_edit_distance=chimeric_stitched_contig_edit_distance,
                                      chimeric_stitched_contig_discordant_reads_cutoff=
                                      chimeric_stitched_contig_discordant_reads_cutoff,
                                      worker_configurer_func=worker_configurer,
                                      counter=counter,
                                      lock=lock,
                                      genes_to_process=genes_to_process,
                                      intronerate=intronerate,
                                      no_padding_supercontigs=no_padding_supercontigs,
                                      keep_intermediate_files=keep_intermediate_files)
                          for gene_name in genes]
        for future in future_results:
            future.add_done_callback(done_callback)

        # As per-gene Exonerate runs complete, read the gene log, log it to the main logger, delete gene log:
        for future in as_completed(future_results):
            try:
                gene_name, prot_length = future.result()
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
            except:  # FIXME make this more specific
                logger.info(f'result is {future.result()}')
                raise

        wait(future_results, return_when="ALL_COMPLETED")

        # Write the 'gene_name', 'prot_length' strings returned by each process to file:
        with open('genes_with_seqs.txt', 'w') as genes_with_seqs_handle:
            for future in future_results:
                try:
                    gene_name, prot_length = future.result()
                    if gene_name and prot_length:
                        genes_with_seqs_handle.write(f'{gene_name}\t{prot_length}\n')
                except ValueError:
                    logger.info(f'result is {future.result()}')
                    raise


def assemble(args):
    """
    Assemble gene, intron, supercontig and paralog sequences via assemble.py

    :param argparse.Namespace args: argparse namespace with subparser options for function assemble()
    :return None: no return value specified; default is None
    """

    # Get a list of read files from args.readfiles (doesn't include any readfile passed in via --unpaired flag):
    readfiles = [os.path.abspath(x) for x in args.readfiles]

    # Check that the --readfiles/-r parameter has been supplied, as it's an expected argument for make_basename():
    if len(readfiles) == 0:
        sys.exit(f'Please provide read files using the parameter --readfiles or -r')

    # Generate a directory for the sample:
    basedir, basename = make_basename(args.readfiles, prefix=args.prefix)

    # Create logger:
    if args.prefix:
        logger = setup_logger(__name__, f'{basename}/{args.prefix}_hybpiper_assemble')
    else:
        logger = setup_logger(__name__, f'{basename}/{os.path.split(readfiles[0])[1].split("_")[0]}_hybpiper_assemble')

    logger.info(f'{"[NOTE]:":10} HybPiper was called with these arguments:\n{" ".join(sys.argv)}\n')

    ####################################################################################################################
    # Check dependencies
    ####################################################################################################################
    if check_dependencies(logger=logger):
        logger.info(f'{"[NOTE]:":10} Everything looks good!')
    else:
        logger.error('ERROR: One or more dependencies not found!')
        return

    ####################################################################################################################
    # Check read and target files
    ####################################################################################################################
    # Check that the target-file and input read files exist and aren't empty:
    for read_file in readfiles:
        if os.path.isfile(read_file) and not os.path.getsize(read_file) == 0:
            logger.debug(f'Input read file {read_file} exists and is not empty, proceeding...')
        else:
            sys.exit(f'Input read file {read_file} does not exist or is empty!')
    if args.unpaired:
        if os.path.isfile(args.unpaired) and not os.path.getsize(args.unpaired) == 0:
            logger.debug(f'Input read file {args.unpaired} exists and is not empty, proceeding...')
        else:
            sys.exit(f'Input read file {args.unpaired} does not exist or is empty!')
    if os.path.isfile(args.targetfile) and not os.path.getsize(args.targetfile) == 0:
        logger.debug(f'Input target/target file {args.targetfile} exists and is not empty, proceeding...')
    else:
        sys.exit(f'Input target/target file {args.targetfile} does not exist or is empty!')

    # If only a single readfile is supplied, set --merged to False regardless of user input:
    if len(readfiles) == 1 and args.merged:
        logger.info(f'{"[NOTE]:":10} The flag --merged has been provided but only a single read file has been '
                    f'supplied. Setting --merged to False.')
        args.merged = False

    # If a file of unpaired reads is provided via the --unpaired parameter and only a single readfile is provided via
    # the -r/--readfiles parameter, exit with an error message:
    if len(readfiles) == 1 and args.unpaired:
        sys.exit(f'{"[ERROR]:":10} You have provided a single file of reads using the -r/--readfiles parameter ('
                 f'{os.path.basename(readfiles[0])}), along with a file of unpaired reads via the --unpaired '
                 f'parameter ({os.path.basename(args.unpaired)}). Please concatenate these two files and provide the '
                 f'single file as input using the -r/--readfiles parameter')

    ####################################################################################################################
    # Read in the target file and read files
    ####################################################################################################################
    if args.targetfile:
        targetfile = os.path.abspath(args.targetfile)
    else:
        print(__doc__)
        return

    # Check that the target file is formatted correctly and translates correctly. If it comprises nucleotides but
    # arg.bwa is false, translate and return path to translated file:
    targetfile = check_targetfile(targetfile, args.bwa, args.targetfile_ambiguity_codes, logger=logger)

    if args.unpaired:
        unpaired_readfile = os.path.abspath(args.unpaired)
    else:
        unpaired_readfile = False
    if len(args.readfiles) < 1:
        logger.error('ERROR: Please specify readfiles with -r')
        return
    if not args.targetfile:
        logger.error('ERROR: Please specify a FASTA file containing target sequences.')
        return

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
    # Map reads to nucleotide targets with BWA
    ####################################################################################################################
    # if args.bwa and not args.no_bwa:
    if args.bwa:
        if args.blast:
            args.blast = False
            if args.unpaired:
                # Note that unpaired_readfile is a single path to the file:
                bwa(unpaired_readfile, targetfile, basename, cpu=args.cpu, unpaired=True, logger=logger)

            # Note that readfiles is a list of one (single-end) or two (paired-end) paths to read files:
            bamfile = bwa(readfiles, targetfile, basename, cpu=args.cpu, logger=logger)
            if not bamfile:
                logger.error(f'{"[ERROR]:":10} Something went wrong with the BWA step, exiting. Check the '
                             f'hybpiper_assemble.log file for sample {basename}!')
                return

            logger.debug(f'bamfile is: {bamfile}')
        else:
            bamfile = f'{basename}.bam'

    ####################################################################################################################
    # Map reads to protein targets with BLASTx
    ####################################################################################################################
    if args.blast:
        if args.unpaired:
            blastx(unpaired_readfile, targetfile, args.evalue, basename, cpu=args.cpu,
                   max_target_seqs=args.max_target_seqs, unpaired=True, logger=logger, diamond=args.diamond,
                   diamond_sensitivity=args.diamond_sensitivity)

        blastx_outputfile = blastx(readfiles, targetfile, args.evalue, basename, cpu=args.cpu,
                                   max_target_seqs=args.max_target_seqs, logger=logger, diamond=args.diamond,
                                   diamond_sensitivity=args.diamond_sensitivity)

        if not blastx_outputfile:
            logger.error(f'{"[ERROR]:":10} Something went wrong with the Blastx step, exiting. Check the '
                         f'hybpiper_assemble.log file for sample {basename}!')
            return
        else:
            blastx_outputfile = f'{basename}.blastx'

    ####################################################################################################################
    # Distribute reads to gene directories from either BLASTx or BWA mapping
    ####################################################################################################################

    if args.distribute:
        pre_existing_fastas = glob.glob('./*/*_interleaved.fasta') + glob.glob('./*/*_unpaired.fasta')
        for fn in pre_existing_fastas:
            os.remove(fn)
        if args.bwa:
            distribute_bwa(bamfile, readfiles, targetfile, target, unpaired_readfile, args.exclude,
                           merged=args.merged, logger=logger)
        else:  # distribute BLASTx results
            distribute_blastx(blastx_outputfile, readfiles, targetfile, target, unpaired_readfile, args.exclude,
                              merged=args.merged, logger=logger)
    if len(readfiles) == 2:
        genes = [x for x in os.listdir('.') if os.path.isfile(os.path.join(x, x + '_interleaved.fasta'))]
    else:
        genes = [x for x in os.listdir('.') if os.path.isfile(os.path.join(x, x + '_unpaired.fasta'))]
    if len(genes) == 0:
        logger.error('ERROR: No genes with BLAST hits! Exiting!')
        return

    ####################################################################################################################
    # Assemble reads using SPAdes
    ####################################################################################################################
    # If the --merged flag is provided, merge reads for SPAdes assembly
    if args.merged:
        logger.info(f'{"[NOTE]:":10} Merging reads for SPAdes assembly')
        for gene in genes:
            interleaved_reads_for_merged = f'{gene}/{gene}_interleaved.fastq'
            logger.debug(f'interleaved_reads_for_merged file is {interleaved_reads_for_merged}\n')
            merged_out = f'{gene}/{gene}_merged.fastq'
            unmerged_out = f'{gene}/{gene}_unmerged.fastq'
            bbmerge_command = f'bbmerge.sh interleaved=true in={interleaved_reads_for_merged} out={merged_out}  ' \
                              f'outu={unmerged_out}'
            try:
                result = subprocess.run(bbmerge_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True)
                logger.debug(f'bbmerge check_returncode() is: {result.check_returncode()}')
                logger.debug(f'bbmerge paired stdout is: {result.stdout}')
                logger.debug(f'bbmerge paired stderr is: {result.stderr}')

            except subprocess.CalledProcessError as exc:
                logger.error(f'bbmerge paired FAILED. Output is: {exc}')
                logger.error(f'bbmerge paired stdout is: {exc.stdout}')
                logger.error(f'bbmerge paired stderr is: {exc.stderr}')
                sys.exit('There was an issue when merging reads. Check read files!')

    if args.assemble:
        if len(readfiles) == 1:
            spades_genelist = spades(genes, cov_cutoff=args.cov_cutoff, cpu=args.cpu, kvals=args.kvals,
                                     paired=False, timeout=args.timeout, logger=logger,
                                     keep_folder=args.keep_intermediate_files, single_cell_mode=args.spades_single_cell)
        elif len(readfiles) == 2:
            if args.merged and not unpaired_readfile:
                spades_genelist = spades(genes, cov_cutoff=args.cov_cutoff, cpu=args.cpu, kvals=args.kvals,
                                         timeout=args.timeout, merged=True, logger=logger,
                                         keep_folder=args.keep_intermediate_files,
                                         single_cell_mode=args.spades_single_cell)
            elif args.merged and unpaired_readfile:
                spades_genelist = spades(genes, cov_cutoff=args.cov_cutoff, cpu=args.cpu, kvals=args.kvals,
                                         timeout=args.timeout, merged=True, unpaired=True, logger=logger,
                                         keep_folder=args.keep_intermediate_files,
                                         single_cell_mode=args.spades_single_cell)
            elif unpaired_readfile and not args.merged:
                spades_genelist = spades(genes, cov_cutoff=args.cov_cutoff, cpu=args.cpu, kvals=args.kvals,
                                         timeout=args.timeout, unpaired=True, logger=logger,
                                         keep_folder=args.keep_intermediate_files,
                                         single_cell_mode=args.spades_single_cell)
            else:
                spades_genelist = spades(genes, cov_cutoff=args.cov_cutoff, cpu=args.cpu, kvals=args.kvals,
                                         timeout=args.timeout, logger=logger, keep_folder=args.keep_intermediate_files,
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

    # return
    if args.exonerate:
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
                                  pool_threads=args.cpu,
                                  logger=logger,
                                  intronerate=args.intronerate,
                                  no_padding_supercontigs=args.no_padding_supercontigs,
                                  keep_intermediate_files=args.keep_intermediate_files)

    ####################################################################################################################
    # Collate all stitched contig and putative chimera read reports
    ####################################################################################################################
    logger.info(f'\n{"[NOTE]:":10} Generated sequences from {len(open("genes_with_seqs.txt").readlines())} genes!')

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
    logger.info(f'{"[NOTE]:":10} WARNING: Potential long paralogs detected for {len(paralog_warnings_long)} genes!')

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
    logger.info(f'{"[NOTE]:":10} WARNING: Potential paralogs detected via contig depth for'
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


def parse_arguments():
    """
    Creates main parser and add subparsers. Parses command line arguments

    :return argparse.Namespace arguments: arguments for the given command/subcommand
    """
    parser = argparse.ArgumentParser(prog='hybpiper', description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='To view parameters and help for a subcommand, use e.g. "assemble '
                                            '--help"')
    parser.add_argument('--version', '-v',
                        dest='version',
                        action='version',
                        version='%(prog)s 2.0rc build 5',
                        help='Print the HybPiper version number.')
    parser.add_argument('--check_dependencies',
                        dest='check_dependencies',
                        action='store_true',
                        default=False,
                        help='Run the check for all pipeline dependencies and exit')

    # Add subparsers:
    subparsers = parser.add_subparsers(title='Subcommands for HybPiper', description='Valid subcommands:')
    parser_assemble = hybpiper_subparsers.add_assemble_parser(subparsers)
    parser_stats = hybpiper_subparsers.add_stats_parser(subparsers)
    parser_retrieve_sequences = hybpiper_subparsers.add_retrieve_sequences_parser(subparsers)
    parser_paralog_retriever = hybpiper_subparsers.add_paralog_retriever_parser(subparsers)
    parser_gene_recovery_heatmap = hybpiper_subparsers.add_gene_recovery_heatmap_parser(subparsers)

    # Set functions for subparsers:
    parser_assemble.set_defaults(func=assemble)
    parser_stats.set_defaults(func=hybpiper_stats_main)
    parser_retrieve_sequences.set_defaults(func=retrieve_sequences_main)
    parser_paralog_retriever.set_defaults(func=paralog_retriever_main)
    parser_gene_recovery_heatmap.set_defaults(func=gene_recovery_heatmap_main)

    # Parse and return all arguments:
    arguments = parser.parse_args()

    # Get the run directory containing the assemble.py module:
    run_dir = os.path.realpath(os.path.split(sys.argv[0])[0])

    # If the flag --check_dependencies was provided, check all HybPiper dependencies and exit:
    if arguments.check_dependencies:
        if check_dependencies():
            sys.exit(f'{"[NOTE]:":10} Everything looks good!')
        else:
            sys.exit(f'\n{"[ERROR]:":10} One or more dependencies not found!')

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

################################################## END OF SCRIPT #######################################################
