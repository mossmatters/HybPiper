#!/usr/bin/env python

"""
HybPiper Version 1.4 release candidate (February 2022)

***NOTES ON VERSION 1.4***

After installation of the pipeline, all pipeline commands are now accessed via the main command 'hybpiper',
followed by a subcommand to run different parts of the pipeline. The available subcommands can be viewed by typing
'hybpiper -h' or 'hybpiper --help'. They are:

    assemble            Assemble gene, intron, and supercontig sequences
    stats               Gather statistics about the HybPiper run(s)
    retrieve_sequences  Retrieve sequences generated from multiple runs of HybPiper
    recovery_heatmap    Create a gene recovery heatmap for the HybPiper run
    paralog_retriever   Retrieve paralog sequences for a given gene, for all samples

To view available parameters and help for any subcommand, simply type e.g. 'hybpiper assemble -h'.

NOTE: the command/script 'read_first.py' no longer exists, and has been replaced by the subcommand 'assemble'. So,
if you had previously run 'reads_first.py' on a sample using the command e.g.:

    python path_to/reads_first.py -b test_targets.fasta -r NZ281_R*_test.fastq --prefix NZ281 --bwa

...this is now replaced by the command:

    hybpiper assemble -b test_targets.fasta -r NZ281_R*_test.fastq --prefix NZ281 --bwa

NOTE: the recovery of introns and supercontigs, previously achieved via the script 'intronerate.py',
is now incorporated in to the 'hybpiper assemble' command. It can be enabled using the flag
'--run_intronerate', e.g.:

    hybpiper assemble -b test_targets.fasta -r NZ281_R*_test.fastq --prefix NZ281 --bwa --run_intronerate

NOTE: the command/script 'get_seq_lengths.py' no longer exists, and this functionality has been incorporated in to
the command 'hybpiper stats'. The sequence length details that were previously printed to screen are now written to
the file 'seq_lengths.tsv'. Similarity, the stats details that were previously written to screen by
'hybpiper_stats.py' are now written to the file 'hybpiper_stats.tsv'.

For full details of all commands and changes, please reads the wiki page at **LINK** and the changelog at **LINK**.
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

# Import non-standard-library modules:
try:
    import Bio
except ImportError:
    sys.exit(f"Required Python package 'Bio' not found. Is it installed for the Python used to run HybPiper?")
try:
    import progressbar
except ImportError:
    sys.exit(f"Required Python package 'progressbar2' not found. Is it installed for the Python used to run HybPiper?")

# Check that user has the minimum required version of Biopython (1.80):
biopython_version_print = pkg_resources.get_distribution('biopython').version
biopython_version = [int(value) for value in re.split('[.]', biopython_version_print)[:2]]
if biopython_version[0:2] < [1, 80]:
    sys.exit(f"HybPiper required Biopython version 1.80 or above. You are using version {biopython_version_print}. "
             f"Please update your Biopython for the Python use to run HybPiper!")

# Import HybPiper modules required for assemble.py:
import distribute_reads_to_targets_bwa
import distribute_reads_to_targets
import distribute_targets
import spades_runner
import exonerate_hits
import hybpiper_stats
import retrieve_sequences
import paralog_retriever
import gene_recovery_heatmap


# f-strings will produce a 'SyntaxError: invalid syntax' error if not supported by Python version:
f'Must be using Python 3.6 or higher.'


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


def py_which(cmd, mode=os.F_OK | os.X_OK, path=None):
    """
    Given a command, mode, and a PATH string, return the path which conforms to the given mode on the PATH,
    or None if there is no such file. `mode` defaults to os.F_OK | os.X_OK. `path` defaults to the result
    of os.environ.get("PATH"), or can be overridden with a custom search path.

    :param str cmd: executable name to search for
    :param int mode: bitwise OR of integers from os.F_OK (existence of path) and os.X_OK (executable) access checks
    :param str path: None, or contents of $PATH variable (as recovered in function body)
    :return None or path of executable
    """

    # Check that a given file can be accessed with the correct mode. Additionally check that `file` is not a
    # directory, as on Windows directories pass the os.access check.
    def _access_check(fn, mode):
        return (os.path.exists(fn) and os.access(fn, mode)
                and not os.path.isdir(fn))

    # If we're given a path with a directory part, look it up directly rather than referring to PATH directories.
    # This includes checking relative to the current directory, e.g. ./script
    if os.path.dirname(cmd):
        if _access_check(cmd, mode):
            return cmd
        return None

    if path is None:
        path = os.environ.get('PATH', os.defpath)
    if not path:
        return None
    path = path.split(os.pathsep)

    if sys.platform == 'win32':
        # The current directory takes precedence on Windows.
        if not os.curdir in path:
            path.insert(0, os.curdir)

        # PATHEXT is necessary to check on Windows.
        pathext = os.environ.get('PATHEXT', '').split(os.pathsep)
        # See if the given file matches any of the expected path extensions. This will allow us to short circuit when
        # given "python.exe".  If it does match, only test that one, otherwise we have to try others.
        if any([cmd.lower().endswith(ext.lower()) for ext in pathext]):
            files = [cmd]
        else:
            files = [cmd + ext for ext in pathext]
    else:
        # On other platforms you don't have things like PATHEXT to tell you what file suffixes are executable,
        # so just pass on cmd as-is.
        files = [cmd]

    seen = set()
    for dir in path:
        normdir = os.path.normcase(dir)
        if normdir not in seen:
            seen.add(normdir)
            for thefile in files:
                name = os.path.join(dir, thefile)
                if _access_check(name, mode):
                    return name
    return None


def check_dependencies(logger=None):
    """
    Checks for the presence of executables and Python packages. Returns a boolean.

    :param logging.Logger logger: a logger object
    return: bool everything_is_awesome: True if all dependencies are found and are executable
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

    logger.info(f'{"[NOTE]:":10} Checking for external dependencies:\n')
    everything_is_awesome = True
    for e in executables:
        e_loc = py_which(e)
        if e_loc:
            logger.info(f'{e:20} found at {e_loc}')
        else:
            logger.info(f'{e:20} not found in your $PATH!')
            everything_is_awesome = False
    logger.info('')
    return everything_is_awesome


def check_baitfile(baitfile, using_bwa, logger=None):
    """
    Checks bait-file fasta header formatting ("taxon*-unique_gene_ID").
    Reports the number of unique genes (each can have multiple representatives) in the baitfile.
    Performs a quick detection of whether baitfile is DNA or amino-acid.
    Checks that seqs in bait file can be translated from the first codon position in the forwards frame (multiple of
    three, no unexpected stop codons), and logs a warning if not.

    :param str baitfile: path to the baitfile
    :param bool using_bwa: True if the --bwa flag is used; a nucleotide target file is expected in this case
    :param logging.Logger logger: a logger object
    :return: None
    """

    # Check bait-file fasta header formatting:
    logger.info(f'{"[NOTE]:":10} Checking baitfile FASTA header formatting...')
    gene_lists = defaultdict(list)
    with open(baitfile, 'r') as bait_file:
        seqs = list(Bio.SeqIO.parse(bait_file, 'fasta'))
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
        logger.error(f'{"[ERROR!]:":10} The following sequences in your baitfile have incorrectly formatted fasta '
                     f'headers:\n')
        fill = textwrap.fill(f'{seq_list}')
        logger.info(textwrap.indent(fill, ' ' * 11))
        logger.info('')
        sys.exit(1)  # baitfile fasta header formatting should be fixed!
    else:
        logger.info(f'{"[NOTE]:":10} The baitfile FASTA header formatting looks good!')

    # Check for duplicated genes:
    duplicated_genes = []
    for gene, gene_count in check_for_duplicate_genes_dict.items():
        if gene_count > 1:
            duplicated_genes.append(gene)
    if duplicated_genes:
        gene_list = ' '.join(duplicated_genes)
        logger.error(f'{"[ERROR!]:":10} The following sequences in your baitfile occur more than once:\n')
        fill = textwrap.fill(f'{gene_list}')
        logger.info(textwrap.indent(fill, ' ' * 11))
        logger.error(f'\nPlease remove duplicate genes before running HybPiper!')
        sys.exit(1)  # duplicate genes in baitfile should be removed!

    # Report the number of unique genes represented in the baitfile:
    logger.info(f'{"[NOTE]:":10} The baitfile contains at least one sequence for {len(gene_lists)} '
                f'unique genes.')

    # Quick detection of whether baitfile is DNA or amino-acid:
    dna = set('ATCGN')
    if os.path.isfile(baitfile):
        with open(baitfile) as bf:
            header = bf.readline()  # skip the first fasta header
            seqline = bf.readline().rstrip().upper()
            if using_bwa and set(seqline) - dna:
                sys.exit(f'{"[ERROR]:":10} Characters other than ACTGN found in first line. You need a nucleotide bait '
                         f'file for BWA!')
            elif not using_bwa and not set(seqline) - dna:
                sys.exit(f'{"ERROR:":10} Only ATCGN characters found in first line. You need a protein bait file for '
                         f'BLASTx!')

    # Check that seqs in bait file can be translated from the first codon position in the forwards frame:
    if using_bwa:
        seqs_needed_padding_dict = defaultdict(list)
        seqs_with_stop_codons_dict = defaultdict(list)

        for seq in seqs:
            gene_name = seq.name.split('-')[-1]
            sequence, needed_padding = distribute_targets.pad_seq(seq)
            translated_seq = sequence.seq.translate()
            num_stop_codons = translated_seq.count('*')

            if needed_padding:
                seqs_needed_padding_dict[gene_name].append(seq)

            if num_stop_codons == 1 and re.search('[*]', str(translated_seq)[-5:]):
                logger.debug(f'Translated sequence {seq.name} contains a single stop codon in the last 5 amino-acids, '
                             f'proceeding...')
            elif num_stop_codons >= 1:
                seqs_with_stop_codons_dict[gene_name].append(seq)

        if seqs_with_stop_codons_dict:
            seq_list = [seq.name for gene_name, bait_file_Sequence_list in seqs_with_stop_codons_dict.items() for seq
                        in bait_file_Sequence_list]
            logger.info(f'{"[WARN!]:":10} There are {len(seq_list)} sequences in your baitfile that contain unexpected '
                        f'stop codons when translated in the first forwards frame. \n{" " * 11}If your baitfile '
                        f'contains only protein-coding sequences, please check these sequences. \n{" " * 11}Sequence '
                        f'names can be found in the sample log file.\n')
            logger.debug(f'Bait/target file sequences with unexpected stop codons: {seq_list}')

        if seqs_needed_padding_dict:
            seq_list = [seq.name for gene_name, bait_file_Sequence_list in seqs_needed_padding_dict.items() for seq
                        in bait_file_Sequence_list]
            logger.info(f'{"[WARN!]:":10} There are {len(seq_list)} sequences in your baitfile that are not multiples '
                        f'of three. \n{" " * 11}If your baitfile contains only protein-coding sequences, please check '
                        f'these sequences. \n{" " * 11}Sequence names can be found in the sample log file.\n')
            logger.debug(f'Bait/target file sequences that are not multiples of three: {seq_list}')


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


def bwa(readfiles, baitfile, basename, cpu, unpaired=False, logger=None):
    """
    Conduct a BWA search of input reads against the baitfile.

    :param list readfiles: one or more read files used as input to the pipeline
    :param str baitfile: path to baitfile (i.e. the target file)
    :param str basename: directory name for sample
    :param int cpu: number of threads/cpus to use for BWA mapping
    :param bool unpaired: True if an unpaired file has been provided, False if not
    :param logging.Logger logger: a logger object
    :return: None, or the path to the *.bam output file from BWA alignment of sample reads to the bait file
    """

    if os.path.isfile(baitfile):
        logger.warning(f'BAITFILE is: baitfile')
        logger.warning(f"os.path.split(baitfile)[0] + '.amb' is: {os.path.split(baitfile)[0] + '.amb'}")
        if os.path.isfile(os.path.split(baitfile)[0] + '.amb'):
            db_file = baitfile
            logger.warning(f'db_file is: {db_file}')
        else:
            logger.info(f'{"[NOTE]:":10} Making nucleotide bwa index in current directory.')
            baitfiledir = os.path.split(baitfile)[0]
            if baitfiledir:
                if os.path.realpath(baitfiledir) != os.path.realpath('.'):
                    shutil.copy(baitfile, '.')
            db_file = os.path.split(baitfile)[1]
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
        logger.error(f'ERROR: Cannot find baitfile at: {baitfile}')
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


def blastx(readfiles, baitfile, evalue, basename, cpu=None, max_target_seqs=10, unpaired=False, logger=None,
           diamond=False, diamond_sensitivity=False):
    """
    Creates a blast database from the full protein target file, and performs BLASTx searches of sample
    nucleotide read files against the protein database.

    :param list/str readfiles: list of paired read files OR  path to unpaired readfile if unpaired is True
    :param str baitfile: path to baitfile (i.e. the target file)
    :param float evalue: evalue to use for BLASTx searches
    :param str basename: directory name for sample
    :param int cpu: number of threads/cpus to use for BLASTx searches
    :param int max_target_seqs: maximum target sequences specified for BLASTx searches
    :param str/bool unpaired: a path if an unpaired file has been provided, boolean False if not
    :param logging.Logger logger: a logger object
    :param bool diamond: if True use DIAMOND instead of BLASTX
    :param bool/str diamond_sensitivity: sensitivity to use for DIAMOND. Default is False; uses default DIAMOND
    :return: None, or path to *.blastx output file from DIAMOND/BLASTx searches of sample reads vs baitfile
    """

    if os.path.isfile(baitfile):
        logger.info(f'{"[NOTE]:":10} Making protein blastdb in current directory.')
        if os.path.split(baitfile)[0]:
            shutil.copy(baitfile, '.')
        db_file = os.path.split(baitfile)[1]
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
        logger.error(f'Cannot find baitfile at: {baitfile}')
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


def distribute_blastx(blastx_outputfile, readfiles, baitfile, target=None, unpaired_readfile=None, exclude=None,
                      merged=False, logger=None):
    """
    When using blastx, distribute sample reads to their corresponding target file hits.

    :param str blastx_outputfile: tabular format output file of BLASTx search
    :param list readfiles: one or more read files used as input to the pipeline
    :param str baitfile: path to baitfile (i.e. the target file)
    :param str target: specific target(s) to use. Tab-delimited file (one gene per line) or single seq name
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
        target_string = f'--target {target}'
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
    distribute_targets.distribute_targets(baitfile, delim='-', besthits=besthits, translate=False, target=target_string)
    return None


def distribute_bwa(bamfile, readfiles, baitfile, target=None, unpaired_readfile=None, exclude=None, merged=False,
                   logger=None):
    """
    When using BWA mapping, distribute sample reads to their corresponding target file gene matches.

    Distribute the 'best' target file sequence (translated if necessary) to each gene directory.

    :param str bamfile: *.bam output file from BWA alignment of sample reads to the bait file
    :param list readfiles: one or more read files used as input to the pipeline
    :param str baitfile: path to baitfile (i.e. the target file)
    :param str target: specific target(s) to use. Tab-delimited file (one gene per line) or single seq name
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
    if target:
        target_string = f'--target {target}'
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

    besthits = distribute_targets.tailored_target_bwa(bamfile, unpaired_bool, exclude_string)
    distribute_targets.distribute_targets(baitfile, delim='-', besthits=besthits, translate=True, target=target_string)
    return None


def spades(genes, cov_cutoff=8, cpu=None, paired=True, kvals=None, timeout=None, unpaired=False,
           merged=False, logger=None, keep_folder=False):
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
    :param bool keep_folder: If True, don't delete the SPAdes assembly folder after contig recovery
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
                                                 merged=merged)
    logger.info(f'{"[NOTE]:":10} Finished running initial SPAdes assemblies for all genes with reads!')
    if len(spades_failed) > 0:
        with open('failed_spades.txt', 'w') as failed_spadefile:
            failed_spadefile.write('\n'.join(spades_failed))

        spades_duds = spades_runner.rerun_spades('failed_spades.txt', cov_cutoff=cov_cutoff, cpu=cpu)
        logger.info(f'{"[NOTE]:":10} Finished re-running SPAdes assemblies for genes with unsuccessful initial '
                    f'assemblies!')

        if len(spades_duds) == 0:
            logger.info(f'{"[NOTE]:":10} All SPAdes re-runs completed successfully!')
        # else:
        #     logger.error(f'{"[WARN!]:":10} SPAdes re-runs failed for genes {" ".join(spades_duds)}')

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
              nosupercontigs=False,
              bbmap_memory=1,
              bbmap_subfilter=7,
              bbmap_threads=2,
              chimeric_supercontig_edit_distance=5,
              chimeric_supercontig_discordant_reads_cutoff=5,
              worker_configurer_func=None,
              counter=None,
              lock=None,
              genes_to_process=0,
              intronerate=False,
              no_padding_supercontigs=False):
    """
    :param str gene_name: name of a gene that had at least one SPAdes contig
    :param str basename: directory name for sample
    :param int thresh: percent identity threshold for stitching together Exonerate results
    :param float paralog_warning_min_length_percentage: min % of a contig vs ref protein length for a paralog warning
    :param int depth_multiplier: accept full-length Exonerate hit if coverage depth <depth_multiplier>x next best hit
    :param bool nosupercontigs: if True, don't create supercontigs and just use longest Exonerate hit
    :param int bbmap_memory: GB memory (RAM ) to use for bbmap.sh
    :param int bbmap_subfilter: ban alignments with more than this many substitutions
    :param int bbmap_threads: number of threads to use for BBmap when searching for chimeric supercontigs
    :param int chimeric_supercontig_edit_distance: min num differences for a read pair to be flagged as discordant
    :param int chimeric_supercontig_discordant_reads_cutoff: min num discordant reads pairs to flag a supercontig as chimeric
    :param function worker_configurer_func: function to configure logging to file
    :param multiprocessing.managers.ValueProxy counter:
    :param multiprocessing.managers.AcquirerProxy lock:
    :param int genes_to_process: total number of genes to be processed via Exonerate
    :param bool intronerate: if True, run intronerate (if supercontig also produced)
    :param bool no_padding_supercontigs: if True, don't pad contig joins in supercontigs with stretches if 10 Ns
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

    # Set whether the chimeric supercontigs test will be performed, and whether a file of interleaved reads is found:
    perform_supercontig_chimera_test, path_to_interleaved_fasta = exonerate_hits.set_supercontig_chimera_test(
        nosupercontigs,
        prefix)

    logger.debug(f'perform_supercontig_chimera_test is: {perform_supercontig_chimera_test}')
    logger.debug(f'path_to_interleaved_fasta is: {path_to_interleaved_fasta}')

    # Read the SPAdes contigs and the 'best' protein reference seq into SeqIO dictionaries:
    try:
        spades_assembly_dict, best_protein_ref_dict = exonerate_hits.parse_spades_and_best_reference(
            f'{gene_name}/{gene_name}_contigs.fasta',
            f'{gene_name}/{gene_name}_baits.fasta',
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
    exonerate_text_output = exonerate_hits.initial_exonerate(f'{gene_name}/{gene_name}_baits.fasta',
                                                             f'{gene_name}/{gene_name}_contigs.fasta',
                                                             prefix)
    if exonerate_text_output:  # i.e. if the initial_exonerate DID produce a result
        exonerate_result = exonerate_hits.parse_exonerate_and_get_supercontig(
            exonerate_text_output,
            query_file=f'{gene_name}/{gene_name}_baits.fasta',
            paralog_warning_min_length_percentage=paralog_warning_min_length_percentage,
            thresh=thresh,
            logger=logger,
            prefix=prefix,
            discordant_cutoff=chimeric_supercontig_discordant_reads_cutoff,
            edit_distance=chimeric_supercontig_edit_distance,
            bbmap_subfilter=bbmap_subfilter,
            bbmap_memory=bbmap_memory,
            bbmap_threads=bbmap_threads,
            interleaved_fasta_file=path_to_interleaved_fasta,
            nosupercontigs=nosupercontigs)

        if intronerate and exonerate_result and exonerate_result.hits_filtered_by_pct_similarity_dict:
            logger.debug(f'exonerate_result.hits_subsumed_hits_removed_overlaps_trimmed_dict for gene {gene_name} is:'
                         f' {exonerate_result.hits_subsumed_hits_removed_overlaps_trimmed_dict}')
            exonerate_hits.intronerate(exonerate_result, spades_assembly_dict, logger=logger,
                                       no_padding_supercontigs=no_padding_supercontigs)
    else:
        exonerate_result = False

    with lock:
        counter.value += 1
        sys.stderr.write(f'\r{"[NOTE]:":10} Finished running Exonerate for gene {gene_name}, {counter.value}'
                         f'/{genes_to_process}')

    if not exonerate_text_output or not exonerate_result or not exonerate_result.supercontig_seqrecord:
        return gene_name, None  # return gene_name to that exonerate_hits.py log can be re-logged to main log file

    return gene_name, len(exonerate_result.supercontig_seqrecord)


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
                              nosupercontigs=False,
                              bbmap_memory=1,
                              bbmap_subfilter=7,
                              bbmap_threads=2,
                              chimeric_supercontig_edit_distance=5,
                              chimeric_supercontig_discordant_reads_cutoff=5,
                              logger=None,
                              intronerate=False,
                              no_padding_supercontigs=False):
    """
    Runs the function exonerate() using multiprocessing.

    :param list genes: list of genes that had successful SPAdes runs
    :param str basename: directory name for sample
    :param int thresh: percent identity threshold for stitching together Exonerate results
    :param float paralog_warning_min_length_percentage: min % of a contig vs ref protein length for a paralog warning
    :param int pool_threads: number of threads/cpus to use for the ProcessPoolExecutor pool
    :param int depth_multiplier: accept full-length Exonerate hit if coverage depth <depth_multiplier>x next best hit
    :param bool nosupercontigs: if True, don't create supercontigs and just use longest Exonerate hit
    :param int bbmap_memory: GB memory (RAM ) to use for bbmap.sh
    :param int bbmap_subfilter: ban alignments with more than this many substitutions
    :param int bbmap_threads: number of threads to use for BBmap when searching for chimeric supercontigs
    :param int chimeric_supercontig_edit_distance: min num differences for a read pair to be flagged as discordant
    :param int chimeric_supercontig_discordant_reads_cutoff: min num discordant reads pairs to flag a supercontig as
    chimeric
    :param logging.Logger logger: a logger object
    :param bool intronerate: if True, intronerate will be run (if a gene is constructed from hits with introns)
    :param bool no_padding_supercontigs: if True, don't pad contig joins in supercontigs with stretches if 10 Ns
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
                                      nosupercontigs=nosupercontigs,
                                      bbmap_memory=bbmap_memory,
                                      bbmap_subfilter=bbmap_subfilter,
                                      bbmap_threads=bbmap_threads,
                                      chimeric_supercontig_edit_distance=chimeric_supercontig_edit_distance,
                                      chimeric_supercontig_discordant_reads_cutoff=
                                      chimeric_supercontig_discordant_reads_cutoff,
                                      worker_configurer_func=worker_configurer,
                                      counter=counter,
                                      lock=lock,
                                      genes_to_process=genes_to_process,
                                      intronerate=intronerate,
                                      no_padding_supercontigs=no_padding_supercontigs)
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
                    os.remove(gene_log_file_to_cat)
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

    run_dir = os.path.realpath(os.path.split(sys.argv[0])[0])

    # Generate a directory for the sample:
    basedir, basename = make_basename(args.readfiles, prefix=args.prefix)

    # Get a list of read files from args.readfiles (doesn't include any readfile passed in via --unpaired flag):
    readfiles = [os.path.abspath(x) for x in args.readfiles]

    # Create logger:
    if args.prefix:
        logger = setup_logger(__name__, f'{basename}/{args.prefix}_hybpiper_assemble')
    else:
        logger = setup_logger(__name__, f'{basename}/{os.path.split(readfiles[0])[1].split("_")[0]}_hybpiper_assemble')

    logger.info(f'{"[NOTE]:":10} HybPiper was called with these arguments:\n{" ".join(sys.argv)}\n')

    # Check that the target/bait-file and input read files exist and aren't empty:
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
    if os.path.isfile(args.baitfile) and not os.path.getsize(args.baitfile) == 0:
        logger.debug(f'Input target/bait file {args.baitfile} exists and is not empty, proceeding...')
    else:
        sys.exit(f'Input target/bait file {args.baitfile} does not exist or is empty!')

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
    # Check dependencies
    ####################################################################################################################
    if check_dependencies(logger=logger):
        other_scripts = ['distribute_reads_to_targets.py', 'distribute_reads_to_targets_bwa.py',
                         'distribute_targets.py', 'exonerate_hits.py', 'spades_runner.py']
        for script in other_scripts:
            if os.path.isfile(os.path.join(run_dir, script)):
                logger.debug(f'Found script {script}, continuing...')
            else:
                logger.error(f'ERROR: Script {script} not found! Please make sure it is in the same directory as '
                             f'this one!')
                return
        logger.info(f'{"[NOTE]:":10} Everything looks good!')
    else:
        logger.error('ERROR: One or more dependencies not found!')
        return

    ####################################################################################################################
    # Read in the target file (called baitfile here) and read files
    ####################################################################################################################
    if args.baitfile:
        baitfile = os.path.abspath(args.baitfile)
    else:
        print(__doc__)
        return

    # Check that the bait file is formatted correctly, translates correctly:
    check_baitfile(baitfile, args.bwa, logger=logger)

    if args.unpaired:
        unpaired_readfile = os.path.abspath(args.unpaired)
    else:
        unpaired_readfile = False
    if len(args.readfiles) < 1:
        logger.error('ERROR: Please specify readfiles with -r')
        return
    if not args.baitfile:
        logger.error('ERROR: Please specify a FASTA file containing target sequences.')
        return

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
                bwa(unpaired_readfile, baitfile, basename, cpu=args.cpu, unpaired=True, logger=logger)

            # Note that readfiles is a list of one (single-end) or two (paired-end) paths to read files:
            bamfile = bwa(readfiles, baitfile, basename, cpu=args.cpu, logger=logger)
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
            blastx(unpaired_readfile, baitfile, args.evalue, basename, cpu=args.cpu,
                   max_target_seqs=args.max_target_seqs, unpaired=True, logger=logger, diamond=args.diamond,
                   diamond_sensitivity=args.diamond_sensitivity)

        blastx_outputfile = blastx(readfiles, baitfile, args.evalue, basename, cpu=args.cpu,
                                   max_target_seqs=args.max_target_seqs, logger=logger, diamond=args.diamond,
                                   diamond_sensitivity=args.diamond_sensitivity)

        if not blastx_outputfile:
            logger.error(f'{"[ERROR]:":10} Something went wrong with the Blastx step, exiting. Check the '
                         f'hybpiper_assemble.log file for sample {basename}!')
            return
        else:
            blastx_outputfile = f'{basename}.blastx'

    ####################################################################################################################
    # Distribute reads to gene directories for either BLASTx or BWA mapping
    ####################################################################################################################
    if args.distribute:
        pre_existing_fastas = glob.glob('./*/*_interleaved.fasta') + glob.glob('./*/*_unpaired.fasta')
        for fn in pre_existing_fastas:
            os.remove(fn)
        if args.bwa:
            distribute_bwa(bamfile, readfiles, baitfile, args.target, unpaired_readfile, args.exclude,
                           merged=args.merged, logger=logger)
        else:  # distribute BLASTx results
            distribute_blastx(blastx_outputfile, readfiles, baitfile, args.target, unpaired_readfile, args.exclude,
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
                                     paired=False, timeout=args.timeout, logger=logger, keep_folder=args.keep_spades)
        elif len(readfiles) == 2:
            if args.merged and not unpaired_readfile:
                spades_genelist = spades(genes, cov_cutoff=args.cov_cutoff, cpu=args.cpu, kvals=args.kvals,
                                         timeout=args.timeout, merged=True, logger=logger, keep_folder=args.keep_spades)
            elif args.merged and unpaired_readfile:
                spades_genelist = spades(genes, cov_cutoff=args.cov_cutoff, cpu=args.cpu, kvals=args.kvals,
                                         timeout=args.timeout, merged=True, unpaired=True, logger=logger,
                                         keep_folder=args.keep_spades)
            elif unpaired_readfile and not args.merged:
                spades_genelist = spades(genes, cov_cutoff=args.cov_cutoff, cpu=args.cpu, kvals=args.kvals,
                                         timeout=args.timeout, unpaired=True, logger=logger,
                                         keep_folder=args.keep_spades)
            else:
                spades_genelist = spades(genes, cov_cutoff=args.cov_cutoff, cpu=args.cpu, kvals=args.kvals,
                                         timeout=args.timeout, logger=logger, keep_folder=args.keep_spades)

        else:
            logger.error('ERROR: Please specify either one (unpaired) or two (paired) read files! Exiting!')
            return
        if not spades_genelist:
            logger.error('ERROR: No genes had assembled contigs! Exiting!')
            return

    ####################################################################################################################
    # Run Exonerate on the assembled SPAdes contigs, and Intronerate if flag --run_intronerate is used:
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
                                  nosupercontigs=args.nosupercontigs,
                                  bbmap_memory=args.bbmap_memory,
                                  bbmap_subfilter=args.bbmap_subfilter,
                                  chimeric_supercontig_edit_distance=args.chimeric_supercontig_edit_distance,
                                  chimeric_supercontig_discordant_reads_cutoff=
                                  args.chimeric_supercontig_discordant_reads_cutoff,
                                  bbmap_threads=args.bbmap_threads,
                                  pool_threads=args.cpu,
                                  logger=logger,
                                  intronerate=args.intronerate,
                                  no_padding_supercontigs=args.no_padding_supercontigs)

    ####################################################################################################################
    # Collate all supercontig and putative chimera read reports
    ####################################################################################################################
    logger.info(f'\n{"[NOTE]:":10} Generated sequences from {len(open("genes_with_seqs.txt").readlines())} genes!')

    # Supercontigs:
    collate_supercontig_reports = [x for x in glob.glob(f'*/{basename}/genes_with_supercontigs.csv')]
    with open(f'{basename}_genes_with_supercontigs.csv', 'w') as genes_with_supercontigs_handle:
        for report_file in collate_supercontig_reports:
            with open(report_file, 'r') as report_handle:
                lines = report_handle.readlines()
                genes_with_supercontigs_handle.write('\n'.join(lines))

    # Putative chimeras:
    collate_putative_chimeras_reports = [x for x in glob.glob(f'*/{basename}/putative_chimeric_supercontigs.csv')]
    with open(f'{basename}_genes_derived_from_putative_chimera_supercontigs.csv',
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
                report_line = warning_handle.readline()
                long_paralogs_handle.write(report_line)
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


def add_assemble_parser(subparsers):
    """
    Parser for the main assembly stage of HybPiper i.e. reads_first.

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_assemble = subparsers.add_parser('assemble', help='Assemble gene, intron, and supercontig sequences')
    parser_assemble.add_argument('--readfiles', '-r', nargs='+',
                                 help='One or more read files to start the pipeline. If exactly two are specified, '
                                      'will assume it is paired Illumina reads.',
                                 default=[], required=True)
    parser_assemble.add_argument('--baitfile', '-b',
                                 help='FASTA file containing bait sequences for each gene. If there are multiple '
                                      'baits for a gene, the id must be of the form: >Taxon-geneName',
                                 default=None, required=True)
    group_1 = parser_assemble.add_mutually_exclusive_group()
    group_1.add_argument('--bwa', dest='bwa', action='store_true',
                         help='Use BWA to search reads for hits to target. Requires BWA and a bait file that is '
                              'nucleotides!', default=False)
    group_1.add_argument('--diamond', dest='diamond', action='store_true', help='Use DIAMOND instead of BLASTx',
                         default=False)
    parser_assemble.add_argument('--diamond_sensitivity', choices=['mid-sensitive', 'sensitive', 'more-sensitive',
                                                                   'very-sensitive', 'ultra-sensitive'],
                                 help='Use the provided sensitivity for DIAMOND searches', default=False)
    parser_assemble.add_argument('--no-blast', dest='blast', action='store_false',
                                 help='Do not run the blast step. Downstream steps will still depend on the '
                                      '*_all.blastx file. \nUseful for re-running assembly/exonerate steps with '
                                      'different options.')
    parser_assemble.add_argument('--no-distribute', dest='distribute', action='store_false',
                                 help='Do not distribute the reads and bait sequences to sub-directories.')
    parser_assemble.add_argument('--no-assemble', dest='assemble', action='store_false',
                                 help='Skip the SPAdes assembly stage.')
    parser_assemble.add_argument('--no-exonerate', dest='exonerate', action='store_false',
                                 help='Do not run the Exonerate step, which assembles full length CDS regions and '
                                      'proteins from each gene')
    parser_assemble.add_argument('--cpu', type=int, default=0,
                                 help='Limit the number of CPUs. Default is to use all cores available.')
    parser_assemble.add_argument('--evalue', type=float, default=1e-10,
                                 help='e-value threshold for blastx hits, default: %(default)s')
    parser_assemble.add_argument('--max_target_seqs', type=int, default=10,
                                 help='Max target seqs to save in blast search, default: %(default)s')
    parser_assemble.add_argument('--cov_cutoff', type=int, default=8,
                                 help='Coverage cutoff for SPAdes. default: %(default)s')
    parser_assemble.add_argument('--kvals', nargs='+',
                                 help='Values of k for SPAdes assemblies. SPAdes needs to be compiled to handle '
                                      'larger k-values! Default auto-detection by SPAdes.', default=None)
    parser_assemble.add_argument('--thresh', type=int,
                                 help='Percent identity threshold for retaining Exonerate hits. Default is 55, '
                                      'but increase this if you are worried about contaminant sequences.', default=55)
    parser_assemble.add_argument('--paralog_min_length_percentage', default=0.75, type=float,
                                 help='Minimum length percentage of a contig Exonerate hit vs reference protein '
                                      'length for a paralog warning and sequence to be generated. Default is %('
                                      'default)s')
    parser_assemble.add_argument('--depth_multiplier',
                                 help='Accept any full-length exonerate hit if it has a coverage depth X times the '
                                      'next best hit. Set to zero to not use depth. Default = 10', default=10, type=int)
    parser_assemble.add_argument('--prefix',
                                 help='Directory name for pipeline output, default is to use the FASTQ file name.',
                                 default=None)
    parser_assemble.add_argument('--timeout',
                                 help='Use GNU Parallel to kill long-running processes if they take longer than X '
                                      'percent of average.', default=0, type=int)
    parser_assemble.add_argument('--target',
                                 help='Use this target to align sequences for each gene. Other targets for that gene '
                                      'will be used only for read sorting. Can be a tab-delimited file (one gene per '
                                      'line) or a single sequence name', default=None)
    parser_assemble.add_argument('--unpaired',
                                 help='Include a single FASTQ file with unpaired reads along with the two paired read '
                                      'files',
                                 default=False)
    parser_assemble.add_argument('--exclude',
                                 help='Do not use any sequence with the specified string as a target sequence for '
                                      'exonerate. The sequence will be used for read sorting.', default=None)
    parser_assemble.add_argument('--nosupercontigs', dest='nosupercontigs', action='store_true',
                                 help='Do not create any supercontigs. The longest single Exonerate hit will be used',
                                 default=False)
    parser_assemble.add_argument('--bbmap_memory', default=1, type=int,
                                 help='GB memory (RAM ) to use for bbmap.sh with exonerate_hits.py. Default is 1')
    parser_assemble.add_argument('--bbmap_subfilter', default=7, type=int,
                                 help='Ban alignments with more than this many substitutions. Default is %(default)s')
    parser_assemble.add_argument('--bbmap_threads', default=2, type=int,
                                 help='Number of threads to use for BBmap when searching for chimeric supercontigs. '
                                      'Default is %(default)s')
    parser_assemble.add_argument('--chimeric_supercontig_edit_distance',
                                 help='Minimum number of differences between one read of a read pair vs the '
                                      'supercontig reference for a read pair to be flagged as discordant', default=5,
                                 type=int)
    parser_assemble.add_argument('--chimeric_supercontig_discordant_reads_cutoff',
                                 help='Minimum number of discordant reads pairs required to flag a supercontig as a '
                                      'potential chimera of contigs from multiple paralogs', default=5, type=int)
    parser_assemble.add_argument('--merged', help='For assembly with both merged and unmerged (interleaved) reads',
                                 action='store_true', default=False)
    parser_assemble.add_argument("--run_intronerate",
                                 help='Run intronerate to recover fasta files for supercontigs with introns (if '
                                      'present), and introns-only', action='store_true', dest='intronerate',
                                 default=False)
    parser_assemble.add_argument("--keep_spades_folder",
                                 help='Keep the SPAdes folder for each gene. Default action is to delete it following '
                                      'contig recovery (dramatically reduces the total files number)',
                                 action='store_true', dest='keep_spades', default=False)
    parser_assemble.add_argument("--no_padding_supercontigs",
                                 help='If Intronerate is run, and a supercontig is created by concatenating multiple '
                                      'SPAdes contigs, do not add 10 "N" characters between contig joins. By default, '
                                      'Ns will be added.', action='store_true', dest='no_padding_supercontigs',
                                 default=False)

    # Set defaults for subparser <parser_assemble>:
    parser_assemble.set_defaults(check_depend=False, blast=True, distribute=True, assemble=True, exonerate=True, )

    # Set function for subparser <parser_assemble>:
    parser_assemble.set_defaults(func=assemble)


def add_stats_parser(subparsers):
    """
    Parser for hybpiper_stats, which now includes running get_seq_lengths

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_stats = subparsers.add_parser('stats', help='Gather statistics about the HybPiper run(s)')
    parser_stats.add_argument('baitfile', help='FASTA file containing bait sequences for each gene. If there are '
                                               'multiple baits for a gene, the id must be of the form: >Taxon-geneName')
    parser_stats.add_argument('namelist', help="Text file with names of HybPiper output directories, one per line")
    # parser_stats.add_argument("sequence_type", help="Sequence type (dna or aa) of the baitfile used")
    parser_stats.add_argument("sequence_type", help="Sequence type (dna, aa or supercontig) to recover lengths for",
                              choices=["dna", "DNA", "aa", "AA", "supercontig", "SUPERCONTIG"])
    # parser_stats.add_argument('--blastx_adjustment', dest="blastx_adjustment", action='store_true',
    #                           help="Adjust stats for when blastx is used i.e. protein references, in cases where "
    #                                "get_seq_lengths() has been run with parameter <dna> rather than <aa>",
    #                           default=False)
    # Set function for subparser <parser_stats>:
    parser_stats.set_defaults(func=hybpiper_stats_main)


def add_retrieve_sequences_parser(subparsers):
    """
    Parser for retrieve_sequences

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_retrieve_sequences = subparsers.add_parser('retrieve_sequences', help='Retrieve sequences generated from '
                                                                                 'multiple runs of HybPiper')
    parser_retrieve_sequences.add_argument('targetfile', help="FASTA File containing target sequences")
    parser_retrieve_sequences.add_argument('sample_names', help="Directory containing Hybpiper output OR a file "
                                                                "containing HybPiper output names, one per line")
    parser_retrieve_sequences.add_argument('sequence_type', help="Type of sequence to extract",
                                           choices=["dna", "aa", "intron", "supercontig"])
    parser_retrieve_sequences.add_argument("--hybpiper_dir", help="Specify directory containing HybPiper output")
    parser_retrieve_sequences.add_argument("--fasta_dir", help="Specify directory for output FASTA files")

    # Set function for subparser <parser_retrieve_sequences>:
    parser_retrieve_sequences.set_defaults(func=retrieve_sequences_main)


def add_paralog_retriever_parser(subparsers):
    """
    Parser for paralog_retriever

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_paralog_retriever = subparsers.add_parser('paralog_retriever', help='Retrieve paralog sequences for a '
                                                                               'given gene, for all samples')
    parser_paralog_retriever.add_argument('namelist', help="Text file containing list of HybPiper output directories, "
                                                           "one per line.")
    parser_paralog_retriever.add_argument('gene', help="Name of gene to extract paralogs")
    parser_paralog_retriever.add_argument('--fasta_dir_all', help='Specify directory for output FASTA files (ALL)',
                                          default='paralogs_all')
    parser_paralog_retriever.add_argument('--fasta_dir_no_chimeras', help='Specify directory for output FASTA files ('
                                                                          'no putative chimeric sequences)',
                                          default='paralogs_no_chimeras')

    # Set function for subparser <paralog_retriever>:
    parser_paralog_retriever.set_defaults(func=paralog_retriever_main)


def add_gene_recovery_heatmap_parser(subparsers):
    """
    Parser for gene_recovery_heatmap

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_gene_recovery_heatmap = subparsers.add_parser('recovery_heatmap', help='Create a gene recovery heatmap for '
                                                                                  'the HybPiper run')
    parser_gene_recovery_heatmap.add_argument('seq_lengths_file',
                                              help="filename for the seq_lengths file (output of 'hybpiper "
                                                   "get_seq_lengths')")
    parser_gene_recovery_heatmap.add_argument('--heatmap_filename',
                                              help='filename for the output heatmap, saved by default as a *.png file. '
                                                   'Defaults to "heatmap"', default='heatmap')
    parser_gene_recovery_heatmap.add_argument('--figure_length', type=int,
                                              help='Length dimension (in inches) for the output heatmap file. '
                                                   'Default is automatically calculated based on the number of '
                                                   'genes', default=None)
    parser_gene_recovery_heatmap.add_argument('--figure_height', type=int,
                                              help='height dimension (in inches) for the output heatmap file. '
                                                   'Default is automatically calculated based on the number of '
                                                   'samples', default=None)
    parser_gene_recovery_heatmap.add_argument('--sample_text_size', type=int,
                                              help='Size (in points) for the sample text labels in the output heatmap '
                                                   'file. Default is automatically calculated based on the '
                                                   'number of samples', default=None)
    parser_gene_recovery_heatmap.add_argument('--gene_text_size', type=int,
                                              help='Size (in points) for the gene text labels in the output heatmap '
                                                   'file. Default is automatically calculated based on the '
                                                   'number of genes', default=None)
    parser_gene_recovery_heatmap.add_argument('--heatmap_filetype', choices=['png', 'pdf', 'eps', 'tiff', 'svg'],
                                              help='File type to save the output heatmap image as. Default is *.png',
                                              default='png')
    parser_gene_recovery_heatmap.add_argument('--heatmap_dpi', type=int,
                                              help='Dot per inch (DPI) for the output heatmap image. Default is 300',
                                              default='300')

    # Set function for subparser <parser_gene_recovery_heatmap>:
    parser_gene_recovery_heatmap.set_defaults(func=gene_recovery_heatmap_main)


def parse_arguments():
    """
    Creates main parser and add subparsers. Parses command line arguments

    :return argparse.Namespace arguments: arguments for the given command/subcommand
    """
    parser = argparse.ArgumentParser(prog='hybpiper', description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='To view parameters and help for a subcommand, use e.g. "assemble '
                                            '--help"')

    subparsers = parser.add_subparsers(title='Subcommands for HybPiper', description='Valid subcommands:')
    add_assemble_parser(subparsers)
    add_stats_parser(subparsers)
    add_retrieve_sequences_parser(subparsers)
    add_gene_recovery_heatmap_parser(subparsers)
    add_paralog_retriever_parser(subparsers)

    # Parse and return all arguments:
    arguments = parser.parse_args()

    return arguments


def main():

    if len(sys.argv) == 1:
        print(__doc__)
        sys.exit(1)

    # Parse arguments for the subcommand used:
    args = parse_arguments()

    # Run the function associated with the subcommand:
    args.func(args)


########################################################################################################################
# Run the script
#######################################################################################################################
if __name__ == '__main__':
    main()

################################################## END OF SCRIPT #######################################################
