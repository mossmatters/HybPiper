#!/usr/bin/env python

"""
########################################################################################################################
####################################################### NOTES  #########################################################
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
    check_targetfile    Check the target file for issues, then exit
    fix_targetfile      Fix the target file, then exit
    filter_by_length    Filter the sequences output by command "hybpiper retrieve_sequences" by length
                        (absolute and relative to mean)

To view available parameters and help for any subcommand, simply type e.g. 'hybpiper assemble -h'.

==> NOTE <==
The script 'reads_first.py' no longer exists, and has been replaced by the subcommand 'assemble'. So,
if you had previously run 'reads_first.py' on a sample using the command e.g.:

    python /<path_to>/reads_first.py -b test_targets.fasta -r NZ281_R*_test.fastq --prefix NZ281 --bwa

...this is now replaced by the command:

    hybpiper assemble -t_dna test_targets.fasta -r NZ281_R*_test.fastq --prefix NZ281 --bwa

==> NOTE <==
The recovery of introns and supercontigs, previously achieved via the script 'intronerate.py',
is now incorporated in to the 'hybpiper assemble' command. It is run by default, but can be disabled using the flag
'--no_intronerate'.

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
import sys
import re
import textwrap
import importlib.metadata
import cProfile

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
biopython_version_print = importlib.metadata.version('biopython')
biopython_version = [int(value) for value in re.split('[.]', biopython_version_print)[:2]]
if biopython_version[0:2] < [1, 80]:
    sys.exit(f'HybPiper requires Biopython version 1.80 or above. You are using version {biopython_version_print}. '
             f'Please update your Biopython for the Python installation used to run HybPiper!')

# Import HybPiper modules:
from hybpiper.version import __version__
from hybpiper import assemble
from hybpiper import hybpiper_stats
from hybpiper import retrieve_sequences
from hybpiper import paralog_retriever
from hybpiper import gene_recovery_heatmap
from hybpiper import fix_targetfile
from hybpiper import filter_by_length
from hybpiper import hybpiper_subparsers
from hybpiper import utils


########################################################################################################################
# Define functions
########################################################################################################################
def assemble_main(args):
    """
    Calls the function main() from module assemble

    :param args: argparse namespace with subparser options for function assemble()
    :return: None: no return value specified; default is None
    """

    assemble.main(args)


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

    utils.check_dependencies()


def check_targetfile(args):
    """
    Performs targetfile checks. Does not translate a DNA file; low-complexity checks are performed on the target file
    as provided.

    :param args: argparse namespace with subparser options for function check_targetfile()
    :return: None: no return value specified; default is None
    """

    if args.targetfile_dna:
        targetfile = args.targetfile_dna
        targetfile_type = 'DNA'
    elif args.targetfile_aa:
        targetfile = args.targetfile_aa
        targetfile_type = 'protein'

    utils.check_targetfile(targetfile=targetfile,
                           targetfile_type=targetfile_type,
                           no_terminal_stop_codons=args.no_terminal_stop_codons,
                           sliding_window_size=args.sliding_window_size,
                           complexity_minimum_threshold=args.complexity_minimum_threshold,
                           run_profiler=args.run_profiler,
                           running_as_subcommand=True)


def fix_targetfile_standalone(args):
    """
    Calls the function main() from module fix_targetfile

    :param args: argparse namespace with subparser options for function fix_targetfile_standalone()
    :return: None: no return value specified; default is None
    """

    fix_targetfile.main(args)


def filter_by_length_main(args):
    """
    Calls the function main() from module filter_by_length

    :param args: argparse namespace with subparser options for function filter_by_length_main()
    :return: None: no return value specified; default is None
    """

    filter_by_length.main(args)


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
                         version=f'hybpiper {__version__}',
                         help='Print the HybPiper version number.')

    # Add subparsers:
    subparsers = parser.add_subparsers(title='Subcommands for HybPiper', metavar='')
    parser_assemble = hybpiper_subparsers.add_assemble_parser(subparsers)
    parser_stats = hybpiper_subparsers.add_stats_parser(subparsers)
    parser_retrieve_sequences = hybpiper_subparsers.add_retrieve_sequences_parser(subparsers)
    parser_paralog_retriever = hybpiper_subparsers.add_paralog_retriever_parser(subparsers)
    parser_gene_recovery_heatmap = hybpiper_subparsers.add_gene_recovery_heatmap_parser(subparsers)
    parser_check_dependencies = hybpiper_subparsers.add_check_dependencies_parser(subparsers)
    parser_check_targetfile = hybpiper_subparsers.add_check_targetfile_parser(subparsers)
    parser_fix_targetfile = hybpiper_subparsers.add_fix_targetfile_parser(subparsers)
    parser_filter_by_length = hybpiper_subparsers.add_filter_by_length_parser(subparsers)

    # Set functions for subparsers:
    parser_assemble.set_defaults(func=assemble_main)
    parser_stats.set_defaults(func=hybpiper_stats_main)
    parser_retrieve_sequences.set_defaults(func=retrieve_sequences_main)
    parser_paralog_retriever.set_defaults(func=paralog_retriever_main)
    parser_gene_recovery_heatmap.set_defaults(func=gene_recovery_heatmap_main)
    parser_check_dependencies.set_defaults(func=check_dependencies)
    parser_check_targetfile.set_defaults(func=check_targetfile)
    parser_fix_targetfile.set_defaults(func=fix_targetfile_standalone)
    parser_filter_by_length.set_defaults(func=filter_by_length_main)

    # Parse and return all arguments:
    arguments = parser.parse_args()

    return arguments


def main():

    title = textwrap.dedent(
        fr"""
                                                         T
                                                            T
                                             C  G
     _    _            _       _____      T        G        A
    | |  | |          | |     |  _  \  A              A  A
    | |__| | __    __ | |___  | |_| |  _   _____   _____   _____
    |  __  | \ \  / / |  _  \ |  ___/ | | |  _  \ |  _  | |  _  \
    | |  | |  \ \/ /  | |_| | | |     | | | |_| | |  __/  | |  --
    |_|  |_|   \  /   |_____/ |_|     |_| |  ___/ |_____| |_|
               / /                        | |
              /_/                         |_|

        """
    )

    print(title)

    if len(sys.argv) == 1:
        print(__doc__)
        sys.exit(0)

    # Parse arguments for the command/subcommand used:
    args = parse_arguments()

    # Run the function associated with the subcommand, with or without cProfile:
    if args.run_profiler:
        profiler = cProfile.Profile()
        profiler.enable()
        args.func(args)
        profiler.disable()
        csv = utils.cprofile_to_csv(profiler)

        with open(f'{sys.argv[1]}_cprofile.csv', 'w+') as cprofile_handle:
            cprofile_handle.write(csv)
    else:
        args.func(args)


########################################################################################################################
# Run the script
#######################################################################################################################
if __name__ == '__main__':


    main()

################################################## END OF SCRIPT #######################################################
