#!/usr/bin/env python

"""
HybPiper Version 1.4 release candidate (September 2021)

This script is a wrapper around several scripts in the HybPiper pipeline.
It can check whether you have the appropriate dependencies available (see --check-depend).
It makes sure that the other scripts needed are in the same directory as this one.
Command line options are passed to the other executables.
Unless --prefix is set, output will be put within a directory named after your read files.

To see parameters and help type:

python reads_first.py -h
"""

import argparse
import os
import sys
import importlib
import shutil
import subprocess
import glob
import logging
import distribute_reads_to_targets_bwa
import distribute_reads_to_targets
import distribute_targets
import spades_runner
import datetime


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

    :param string name: name for the logger instance
    :param string log_file: filename for log file
    :param string console_level: logger level for logging to console
    :param string file_level: logger level for logging to file
    :param string logger_object_level: logger level for logger object
    :return: a logger object
    """

    # Get date and time string for log filename:
    date_and_time = datetime.datetime.now().strftime("%Y-%m-%d-%H_%M_%S")

    # Log to file:
    file_handler = logging.FileHandler(f'{log_file}_{date_and_time}.log', mode='w')
    file_handler.setLevel(file_level)
    file_format = logging.Formatter('%(asctime)s - %(filename)s - %(name)s - %(funcName)s - %(levelname)s - %('
                                    'message)s')
    file_handler.setFormatter(file_format)

    # Log to Terminal (stdout):
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
        if not normdir in seen:
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

    python_packages = ['Bio']

    everything_is_awesome = True
    for e in executables:
        e_loc = py_which(e)
        if e_loc:
            logger.info(f'{e} found at {e_loc}')
        else:
            logger.info(f'{e} not found in your $PATH!')
            everything_is_awesome = False

    for p in python_packages:
        try:
            i = importlib.import_module(p)
            logger.info(f'Package {p} successfully loaded!')
        except ImportError:
            logger.error(f'Package {p} not found!')
            everything_is_awesome = False
    return everything_is_awesome


def make_basename(readfiles, prefix=None):
    """
    Unless prefix is set, generate a directory based off the readfiles, using everything up to the first underscore.
    If prefix is set, generate the directory "prefix" and set basename to be the last component of the path.

    :param list readfiles: one or more read files to start the pipeline
    :param str prefix: directory name for pipeline output
    :return: parent directory, directory name
    """

    if prefix:
        if not os.path.exists(prefix):
            os.makedirs(prefix)
        prefixParentDir, prefix = os.path.split(prefix)
        if not prefix:
            # if prefix has a trailing /, prefixParentDir will have the / stripped and prefix will be empty,
            # so try again
            prefix = os.path.split(prefixParentDir)[1]
        return prefixParentDir, prefix

    # --prefix is not set on cmd line;  Write output to subdir in "."
    basename = os.path.split(readfiles[0])[1].split('_')[0]
    if not os.path.exists(basename):
        os.makedirs(basename)
    return '.', basename


def bwa(readfiles, baitfile, basename, cpu, unpaired=False, logger=None):
    """
    Conduct BWA search of reads against the baitfile. Returns an error if the second line of the baitfile contains
    characters other than ACTGN.

    :param list readfiles: one or more read files to start the pipeline
    :param str baitfile: path to baitfile (target file)
    :param str basename: directory name for sample
    :param int cpu: number of threads/cpus to use for BWA mapping
    :param str/bool unpaired: a path if an unpaired file has been provided, False if not
    :param logging.Logger logger: a logger object
    :return: None, or the *.bam output file from BWA alignment of sample reads to the bait file
    """

    dna = set('ATCGN')
    if os.path.isfile(baitfile):
        # Quick detection of whether baitfile is DNA.
        with open(baitfile) as bf:
            header = bf.readline()
            seqline = bf.readline().rstrip().upper()
            if set(seqline) - dna:
                logger.error('ERROR: characters other than ACTGN found in first line. You need a nucleotide bait file '
                             'for BWA!')
                return None

        if os.path.isfile(os.path.split(baitfile)[0] + '.amb'):
            db_file = baitfile
        else:
            logger.info(f'{"[NOTE]:":10} Making nucleotide bwa index in current directory.')
            baitfileDir = os.path.split(baitfile)[0]
            if baitfileDir:
                if os.path.realpath(baitfileDir) != os.path.realpath('.'):
                    shutil.copy(baitfile, '.')
            db_file = os.path.split(baitfile)[1]
            make_bwa_index_cmd = f'bwa index {db_file}'
            logger.info(f'[CMD]: {make_bwa_index_cmd}')

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
        cpu = multiprocessing.cpu_count()

    if len(readfiles) < 3:
        bwa_fastq = ' '.join(readfiles)
    else:
        bwa_fastq = readfiles

    bwa_commands = ['time bwa mem', '-t', str(cpu), db_file, bwa_fastq, ' | samtools view -h -b -S - > ']
    if unpaired:
        bwa_commands.append(f'{basename}_unpaired.bam')
    else:
        bwa_commands.append(f'{basename}.bam')
    full_command = ' '.join(bwa_commands)
    logger.info(f'[CMD]: {full_command}')

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

    return f'{basename}.bam'


def blastx(readfiles, baitfile, evalue, basename, cpu=None, max_target_seqs=10, unpaired=False, logger=None,
           diamond=False, diamond_sensitivity=False):
    """
    Creates a blast database from the complete protein target file, and performs BLASTx searches of sample
    nucleotide read files against the protein database.

    :param list readfiles: one or more read files to start the pipeline
    :param str baitfile: path to baitfile (target file)
    :param float evalue: evalue to use for BLASTx searches
    :param str basename: directory name for sample
    :param int cpu: number of threads/cpus to use for BLASTx searches
    :param int max_target_seqs: maximum target sequences specified for BLASTx searches
    :param str/bool unpaired: a path if an unpaired file has been provided, False if not
    :param logging.Logger logger: a logger object
    :param bool diamond: if True use DIAMOND instead of BLASTX
    :param bool/str diamond_sensitivity: sensitivity to use for DIAMOND. Default is False; uses default DIAMOND
    :return: None, or the *.blastx output file from DIAMOND/BLASTx searches of sample reads against the bait file
    """

    dna = set('ATCGN')
    if os.path.isfile(baitfile):
        # Quick detection of whether baitfile is DNA.
        with open(baitfile) as bf:
            header = bf.readline()
            seqline = bf.readline().rstrip().upper()
            if not set(seqline) - dna:
                logger.info('ERROR: only ATCGN characters found in first line. You need a protein bait file for '
                            'BLASTx!')
                return None

        # if os.path.isfile(os.path.split(baitfile)[0] + '.psq'):
        #     db_file = baitfile
        # else:
        logger.info('Making protein blastdb in current directory.')
        if os.path.split(baitfile)[0]:
            shutil.copy(baitfile, '.')
        db_file = os.path.split(baitfile)[1]
        if diamond:
            logger.info(f'Using DIAMOND instead of BLASTx!')
            if diamond_sensitivity:
                logger.info(f'Using DIAMOND sensitivity "{diamond_sensitivity}"')
            makeblastdb_cmd = f'diamond makedb --in {db_file} --db {db_file}'
        else:
            makeblastdb_cmd = f'makeblastdb -dbtype prot -in {db_file}'
        logger.info(f'[CMD]: {makeblastdb_cmd}')
        try:
            result = subprocess.run(makeblastdb_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True)
            logger.debug(f'makeblastdb/makedb check_returncode() is: {result.check_returncode()}')
            logger.debug(f'makeblastdb/makedb stdout is: {result.stdout}')
            logger.debug(f'makeblastdb/makedb stderr is: {result.stderr}')
        except subprocess.CalledProcessError as exc:
            logger.error(f'makeblastdb/makedb FAILED. Output is: {exc}')
            logger.error(f'makeblastdb/makedb stdout is: {exc.stdout}')
            logger.error(f'makeblastdb/makedb stderr is: {exc.stderr}')
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
        logger.info(f'[CMD]: {full_command}')

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
                pipe_cmd = f"cat {read_file} | \
                awk '{{if(NR % 4 == 1 || NR % 4 == 2) {{sub(/@/, \">\"); print; }} }}'"

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

            logger.info(f'[CMD]: {full_command}')

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

    :param str blastx_outputfile:
    :param list readfiles: one or more read files to start the pipeline
    :param str baitfile: path to baitfile (target file)
    :param str target: specific target(s) to use. Tab-delimited file (one gene per line) or single seq name
    :param str/bool unpaired_readfile: a path if an unpaired file has been provided, False if not
    :param str exclude: specify sequence not to be used as a target sequence for Exonerate
    :param bool merged: if True, write and distribute fastq files for merging with BBmerge.sh (in addition to fasta)
    :param logging.Logger logger: a logger object
    :return: None
    """

    logger.info(f'READFILES are: {readfiles}')

    # Distribute reads to gene directories:
    read_hit_dict_paired = distribute_reads_to_targets.read_sorting(blastx_outputfile)
    logger.info(f'[NOTE]: Unique reads with hits: {len(read_hit_dict_paired)}')
    distribute_reads_to_targets.distribute_reads(readfiles, read_hit_dict_paired, merged=merged)

    if unpaired_readfile:
        blastx_outputfile = blastx_outputfile.replace('.blastx', '_unpaired.blastx')
        read_hit_dict_unpaired = distribute_reads_to_targets.read_sorting(blastx_outputfile)
        distribute_reads_to_targets.distribute_reads([unpaired_readfile], read_hit_dict_unpaired)

    # Distribute the 'best' target file sequence (translated if necessary) to each gene directory:
    if target:
        target_string = f'--target {target}'
    else:
        target_string = None
    # if unpaired:  # CJJ blastX doesn't use unpaired file to calculate best target seq?
    #     unpaired_bool = True
    # else:
    #     unpaired_bool = False
    if exclude:
        exclude_string = f'--exclude {exclude}'
    else:
        exclude_string = None

    besthits = distribute_targets.tailored_target_blast(blastx_outputfile, exclude_string)
    distribute_targets.distribute_targets(baitfile, dirs=True, delim='-', besthits=besthits, translate=False,
                                          target=target_string)
    return None


def distribute_bwa(bamfile, readfiles, baitfile, target=None, unpaired_readfile=None, exclude=None, merged=False,
                   logger=None):
    """
    When using BWA mapping, distribute sample reads to their corresponding target file gene matches.

    Distribute the 'best' target file sequence (translated if necessary) to each gene directory.

    :param str bamfile: *.bam output file from BWA alignment of sample reads to the bait file
    :param list readfiles: one or more read files to start the pipeline
    :param str baitfile: path to baitfile (target file)
    :param str target: specific target(s) to use. Tab-delimited file (one gene per line) or single seq name
    :param str/bool unpaired_readfile: a path if an unpaired file has been provided, False if not
    :param str exclude: specify sequence not to be used as a target sequence for Exonerate
    :param bool merged: if True, write and distribute fastq files for merging with BBmerge.sh (in addition to fasta)
    :param logging.Logger logger: a logger object
    :return: None
    """

    # Distribute reads to gene directories:
    read_hit_dict_paired = distribute_reads_to_targets_bwa.read_sorting(bamfile)
    logger.info(f'[NOTE]: Unique reads with hits: {len(read_hit_dict_paired)}')
    distribute_reads_to_targets_bwa.distribute_reads(readfiles, read_hit_dict_paired, merged=merged)

    if unpaired_readfile:
        up_bamfile = bamfile.replace('.bam', '_unpaired.bam')
        read_hit_dict_unpaired = distribute_reads_to_targets_bwa.read_sorting(up_bamfile)
        distribute_reads_to_targets_bwa.distribute_reads([unpaired_readfile], read_hit_dict_unpaired)

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
    distribute_targets.distribute_targets(baitfile, dirs=True, delim='-', besthits=besthits, translate=True,
                                          target=target_string)
    return None


def spades(genes, cov_cutoff=8, cpu=None, paired=True, kvals=None, timeout=None, unpaired=False,
           merged=False, logger=None):
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
    :return: list spades_genelist: a list of gene names that had successful SPAdes assemblies
    """

    with open('spades_genelist.txt', 'w') as spadesfile:
        spadesfile.write('\n'.join(genes) + '\n')

    if os.path.isfile('spades.log'):
        os.remove('spades.log')
    if os.path.isfile('spades_redo.log'):
        os.remove('spades_redo.log')

    logger.debug(f'args.merged is: {merged}')
    logger.debug(f'args.unpaired is:  {unpaired}')

    if unpaired:  # Create empty unpaired file if it doesn't exist
        for gene in open('spades_genelist.txt'):
            gene = gene.rstrip()
            if os.path.isfile(f'{gene}/{gene}_interleaved.fasta'):
                if not os.path.isfile(f'{gene}/{gene}_unpaired.fasta'):
                    open(f'{gene}/{gene}_unpaired.fasta', 'a').close()

    spades_failed = spades_runner.spades_initial('spades_genelist.txt', cov_cutoff=cov_cutoff, cpu=cpu,
                                                 kvals=kvals, paired=paired, timeout=timeout, unpaired=unpaired,
                                                 merged=merged)
    if len(spades_failed) > 0:
        with open('failed_spades.txt', 'w') as failed_spadefile:
            failed_spadefile.write('\n'.join(spades_failed))

        spades_duds = spades_runner.rerun_spades('failed_spades.txt', cov_cutoff=cov_cutoff, cpu=cpu)

        if len(spades_duds) == 0:
            logger.info('All redos completed successfully!\n')
        else:
            logger.error(f'SPAdes redos failed for genes {" ".join(spades_duds)}')
            sys.exit(1)

    if os.path.isfile('spades_duds.txt'):  # Written by spades_runner.rerun_spades()
        spades_duds = [x.rstrip() for x in open('spades_duds.txt')]
    else:
        spades_duds = []

    spades_genelist = []
    for gene in genes:
        if gene not in set(spades_duds):
            spades_genelist.append(gene)

    with open('exonerate_genelist.txt', 'w') as genefile:
        genefile.write('\n'.join(spades_genelist) + '\n')
    return spades_genelist


def exonerate(genes, basename, run_dir, replace=True, cpu=None, thresh=55, depth_multiplier=0,
              length_pct=100, timeout=None, nosupercontigs=False, memory=1, discordant_reads_edit_distance=7,
              discordant_reads_cutoff=100, paralog_warning_min_cutoff=0.75, bbmap_subfilter=7, logger=None):
    """
    Runs the `exonerate_hits.py.mossmatters script via GNU parallel.

    :param list genes: list of genes that had successful SPAdes runs
    :param str basename: directory name for sample
    :param run_dir: CJJ will be removed after refactor
    :param bool replace: CJJ hardcoded as True, remove after refactor
    :param int cpu: number of threads/cpus to use for GNU Parallel
    :param int thresh: percent identity threshold for stitching together Exonerate hits
    :param int depth_multiplier: if Exonerate hit if it has coverage depth X times the next best hit, accept
    :param int length_pct: Include Exonerate hit if >= long as X percentage of the reference protein length
    :param int timeout: value for GNU parallel --timeout percentage
    :param bool nosupercontigs: If True, no not produce supercontigs; use longest Exonerate hit only
    :param int memory: GB memory (RAM ) to use for bbmap.sh
    :param int discordant_reads_edit_distance: minimum number of edits for discordant read pair flagging
    :param int discordant_reads_cutoff: number of discordant reads for supercontig to be flagged as a chimera
    :param float paralog_warning_min_cutoff: min length % of a contig vs reference protein length for paralog warning
    :param int bbmap_subfilter: ban alignments with more than this many substitutions
    :param logging.Logger logger: a logger object
    :return: None
    """

    if replace:
        for g in genes:
            if os.path.isdir(os.path.join(g, basename)):
                shutil.rmtree(os.path.join(g, basename))
    if len(genes) == 0:
        logger.error(f'ERROR: No genes recovered for {basename}!')
        return 1

    if os.path.isfile('genes_with_seqs.txt'):
        os.remove('genes_with_seqs.txt')

    file_stem = "contigs.fasta"

    parallel_cmd_list = ['time parallel', '--eta', '--joblog parallel.log']
    if cpu:
        parallel_cmd_list.append(f'-j {cpu}')
    if timeout:
        parallel_cmd_list.append(f'--timeout {timeout}%')

    if nosupercontigs:
        logger.info(f'Running Exonerate to generate sequences for {len(genes)} genes, without creating supercontigs')
        exonerate_cmd_list = ['python',
                              f'{run_dir}/exonerate_hits.py.mossmatters',
                              '{}/{}_baits.fasta',
                              '{{}}/{{}}_{}'.format(file_stem),
                              '--prefix {{}}/{}'.format(basename),
                              f'-t {thresh}',
                              f'--depth_multiplier {depth_multiplier}',
                              f'--length_pct {length_pct}',
                              '--nosupercontigs',
                              f'--paralog_warning_min_cutoff {paralog_warning_min_cutoff}',
                              f'--bbmap_subfilter {bbmap_subfilter}',
                              '::::',
                              'exonerate_genelist.txt',
                              '> genes_with_seqs.txt']
    else:
        logger.info(f'Running Exonerate to generate sequences for {len(genes)} genes')
        exonerate_cmd_list = ['python',
                              f'{run_dir}/exonerate_hits.py.mossmatters',
                              '{}/{}_baits.fasta',
                              '{{}}/{{}}_{}'.format(file_stem),
                              '--prefix {{}}/{}'.format(basename),
                              f'-t {thresh}',
                              f'--depth_multiplier {depth_multiplier}',
                              f'--length_pct {length_pct}',
                              f'--memory {memory}',
                              f'--discordant_reads_edit_distance {discordant_reads_edit_distance}',
                              f'--discordant_reads_cutoff {discordant_reads_cutoff}',
                              f'--paralog_warning_min_cutoff {paralog_warning_min_cutoff}',
                              f'--bbmap_subfilter {bbmap_subfilter}',
                              '::::',
                              'exonerate_genelist.txt',
                              "> genes_with_seqs.txt"]

    exonerate_cmd = ' '.join(parallel_cmd_list) + ' ' + ' '.join(exonerate_cmd_list)
    logger.info(exonerate_cmd)
    exitcode = subprocess.call(exonerate_cmd, shell=True)

    if exitcode:
        logger.info(f'exitcode is: {exitcode}')
        logger.error('ERROR: Something went wrong with Exonerate!')
        return exitcode
    return


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--check-depend', dest='check_depend',
                        help='Check for dependencies (executables and Python packages) and exit. May not work at all '
                             'on Windows.', action='store_true')
    group_1 = parser.add_mutually_exclusive_group()
    group_1.add_argument('--bwa', dest='bwa', action='store_true',
                       help='Use BWA to search reads for hits to target. Requires BWA and a bait file that is '
                            'nucleotides!', default=False)
    group_1.add_argument('--diamond', dest='diamond', action='store_true', help='Use DIAMOND instead of BLASTx',
                         default=False)
    parser.add_argument('--diamond_sensitivity', choices=['mid-sensitive', 'sensitive', 'more-sensitive',
                                                          'very-sensitive', 'ultra-sensitive'],
                        help='Use the provided sensitivity for DIAMOND searches', default=False)
    parser.add_argument('--no-blast', dest='blast', action='store_false',
                        help='Do not run the blast step. Downstream steps will still depend on the *_all.blastx file. '
                             '\nUseful for re-runnning assembly/exonerate steps with different options.')
    parser.add_argument('--no-distribute', dest='distribute', action='store_false',
                        help='Do not distribute the reads and bait sequences to sub-directories.')
    parser.add_argument('--no-exonerate', dest='exonerate', action='store_false',
                        help='Do not run the Exonerate step, which assembles full length CDS regions and proteins from '
                             'each gene')
    parser.add_argument('--no-assemble', dest='assemble', action='store_false', help='Skip the SPAdes assembly stage.')
    parser.add_argument('-r', '--readfiles', nargs='+',
                        help='One or more read files to start the pipeline. If exactly two are specified, will assume '
                             'it is paired Illumina reads.',
                        default=[])
    parser.add_argument('-b', '--baitfile',
                        help='FASTA file containing bait sequences for each gene. If there are multiple baits for a '
                             'gene, the id must be of the form: >Taxon-geneName',
                        default=None)
    parser.add_argument('--cpu', type=int, default=0,
                        help='Limit the number of CPUs. Default is to use all cores available.')
    parser.add_argument('--evalue', type=float, default=1e-10,
                        help='e-value threshold for blastx hits, default: %(default)s')
    parser.add_argument('--max_target_seqs', type=int, default=10,
                        help='Max target seqs to save in blast search, default: %(default)s')
    parser.add_argument('--cov_cutoff', type=int, default=8, help='Coverage cutoff for SPAdes. default: %(default)s')
    parser.add_argument('--kvals', nargs='+',
                        help='Values of k for SPAdes assemblies. SPAdes needs to be compiled to handle larger k-values!'
                             ' Default auto-detection by SPAdes.', default=None)
    # CJJ: changed from 65 to 55 as I noticed cases with real hits falling beneath cutoff threshold:
    parser.add_argument('--thresh', type=int,
                        help='Percent Identity Threshold for stitching together exonerate results. Default is 55, but '
                             'increase this if you are worried about contaminant sequences.', default=55)

    parser.add_argument('--paralog_warning_min_length_percentage', default=0.75, type=float,
                        help='Minimum length percentage of a contig vs reference protein length for a paralog warning '
                             'to be generated. Default is %(default)s')
    parser.add_argument('--length_pct', help='Include an exonerate hit if it is at least as long as X percentage of '
                                             'the reference protein length. Default = 90%%', default=90, type=int)
    parser.add_argument('--depth_multiplier',
                        help='Accept any full-length exonerate hit if it has a coverage depth X times the next best '
                             'hit. Set to zero to not use depth. Default = 10', default=10, type=int)
    parser.add_argument('--prefix', help='Directory name for pipeline output, default is to use the FASTQ file name.',
                        default=None)
    parser.add_argument('--timeout',
                        help='Use GNU Parallel to kill long-running processes if they take longer than X percent of '
                             'average.', default=0, type=int)
    parser.add_argument('--target',
                        help='Use this target to align sequences for each gene. Other targets for that gene will be '
                             'used only for read sorting. Can be a tab-delimited file (one gene per line) or a single '
                             'sequence name', default=None)
    parser.add_argument('--unpaired',
                        help='Include a single FASTQ file with unpaired reads along with the two paired read files',
                        default=False)
    parser.add_argument('--exclude',
                        help='Do not use any sequence with the specified string as a target sequence for exonerate. '
                             'The sequence will be used for read sorting.', default=None)
    parser.add_argument('--nosupercontigs', dest='nosupercontigs', action='store_true',
                        help='Do not create any supercontigs. The longest single Exonerate hit will be used',
                        default=False)
    parser.add_argument('--memory', help='GB memory (RAM ) to use for bbmap.sh with exonerate_hits.py.mossmatters. Default is 1',
                        default=1, type=int)
    parser.add_argument('--bbmap_subfilter', default=7, type=int,
                        help='Ban alignments with more than this many substitutions. Default is %(default)s')
    parser.add_argument('--discordant_reads_edit_distance',
                        help='Minimum number of differences between one read of a read pair vs the supercontig '
                             'reference for a read pair to be flagged as discordant', default=5, type=int)
    parser.add_argument('--discordant_reads_cutoff',
                        help='minimum number of discordant reads pairs required to flag a supercontigs as a potential '
                             'hybrid of contigs from multiple paralogs', default=5, type=int)
    parser.add_argument('--merged', help='For assembly with both merged and unmerged (interleaved) reads',
                        action='store_true', default=False)

    parser.set_defaults(check_depend=False, blast=True, distribute=True, assemble=True, exonerate=True, )

    if len(sys.argv) == 1:
        print(__doc__)
        sys.exit(1)

    args = parser.parse_args()
    run_dir = os.path.realpath(os.path.split(sys.argv[0])[0])

    # Generate a directory for the sample:
    basedir, basename = make_basename(args.readfiles, prefix=args.prefix)

    # Get a list of read files from args.readfiles (doesn't include any readfile passed in via --unpaired flag):
    readfiles = [os.path.abspath(x) for x in args.readfiles]

    # Create logger:
    if args.prefix:
        logger = setup_logger(__name__, f'{basename}/{args.prefix}_reads_first')
    else:
        logger = setup_logger(__name__, f'{basename}/{os.path.split(readfiles[0])[1].split("_")[0]}_reads_first')

    logger.info(f'HybPiper was called with these arguments:\n{" ".join(sys.argv)}\n')

    # If only a single readfile is supplied, set --merged to False regardless of user input:
    if len(readfiles) == 1 and args.merged:
        logger.info(f'The flag --merged has been provided but only a single read file has been supplied. Setting '
                    f'--merged to False.')
        args.merged = False

    ####################################################################################################################
    # Check dependencies
    ####################################################################################################################
    if args.check_depend:
        if check_dependencies(logger=logger):
            other_scripts = ['distribute_reads_to_targets.py', 'distribute_reads_to_targets_bwa.py',
                             'distribute_targets.py', 'exonerate_hits.py.mossmatters', 'spades_runner.py']
            for script in other_scripts:
                if os.path.isfile(os.path.join(run_dir, script)):
                    logger.debug(f'Found script {script}, continuing...')
                else:
                    logger.error(f'ERROR: Script {script} not found! Please make sure it is in the same directory as '
                                 f'this one!')
                    return
            logger.info('Everything looks good!')
            return
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
    if args.bwa:
        if args.blast:
            args.blast = False
            if args.unpaired:
                bwa(unpaired_readfile, baitfile, basename, cpu=args.cpu, unpaired=True, logger=logger)

            bamfile = bwa(readfiles, baitfile, basename, cpu=args.cpu, logger=logger)
            if not bamfile:
                logger.error(f'ERROR: Something went wrong with the BWA step, exiting. Check the reads_first.log '
                             f'file for sample {basename}!')
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
            logger.error(f'ERROR: Something went wrong with the Blastx step, exiting. Check the reads_first.log file '
                         f'for sample {basename}!')
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
        logger.info(f'Merging reads for SPAdes assembly')
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
                                     paired=False, timeout=args.timeout, logger=logger)
        elif len(readfiles) == 2:
            if args.merged and not unpaired_readfile:
                spades_genelist = spades(genes, cov_cutoff=args.cov_cutoff, cpu=args.cpu, kvals=args.kvals,
                                         timeout=args.timeout, merged=True, logger=logger)
            elif args.merged and unpaired_readfile:
                spades_genelist = spades(genes, cov_cutoff=args.cov_cutoff, cpu=args.cpu, kvals=args.kvals,
                                         timeout=args.timeout, merged=True, unpaired=True, logger=logger)
            elif unpaired_readfile and not args.merged:
                spades_genelist = spades(genes, cov_cutoff=args.cov_cutoff, cpu=args.cpu, kvals=args.kvals,
                                         timeout=args.timeout, unpaired=True, logger=logger)
            else:
                spades_genelist = spades(genes, cov_cutoff=args.cov_cutoff, cpu=args.cpu, kvals=args.kvals,
                                         timeout=args.timeout, logger=logger)

        else:
            logger.error('ERROR: Please specify either one (unpaired) or two (paired) read files! Exiting!')
            return
        if not spades_genelist:
            logger.error('ERROR: No genes had assembled contigs! Exiting!')
            return

    ####################################################################################################################
    # Run Exonerate on the assembled SPAdes contigs
    ####################################################################################################################
    if args.exonerate:
        genes = [x.rstrip() for x in open('exonerate_genelist.txt').readlines()]
        exitcode = exonerate(genes,
                             basename,
                             run_dir,
                             cpu=args.cpu,
                             thresh=args.thresh,
                             length_pct=args.length_pct,
                             depth_multiplier=args.depth_multiplier,
                             timeout=args.timeout,
                             nosupercontigs=args.nosupercontigs,
                             memory=args.memory,
                             discordant_reads_edit_distance=args.discordant_reads_edit_distance,
                             discordant_reads_cutoff=args.discordant_reads_cutoff,
                             paralog_warning_min_cutoff=args.paralog_warning_min_length_percentage,
                             bbmap_subfilter=args.bbmap_subfilter,
                             logger=logger)
        if exitcode:
            return

    ####################################################################################################################
    # Collate all supercontig and discordant read reports into one file
    ####################################################################################################################
    collate_supercontig_reports = f'find .  -name "genes_with_supercontigs.csv" -exec cat {{}} \; | tee ' \
                                  f'{basename}_genes_with_supercontigs.csv'
    subprocess.call(collate_supercontig_reports, shell=True)
    collate_discordant_supercontig_reports = f'find .  -name "supercontigs_with_discordant_readpairs.csv" ' \
                                             f'-exec cat {{}} \; | tee ' \
                                             f'{basename}_supercontigs_with_discordant_reads.csv'
    subprocess.call(collate_discordant_supercontig_reports, shell=True)

    ####################################################################################################################
    # Report paralog warning and write a paralog warning file
    ####################################################################################################################
    logger.info(f'Generated sequences from {len(open("genes_with_seqs.txt").readlines())} genes!')
    paralog_warnings = [x for x in os.listdir('.') if os.path.isfile(os.path.join(x, basename, 'paralog_warning.txt'))]
    with open('genes_with_paralog_warnings.txt', 'w') as pw:
        pw.write('\n'.join(paralog_warnings))
    logger.info(f'WARNING: Potential paralogs detected for {len(paralog_warnings)} genes!')


########################################################################################################################
# Run the script
#######################################################################################################################
if __name__ == '__main__':
    main()

################################################## END OF SCRIPT #######################################################
