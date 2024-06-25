#!/usr/bin/env python

"""
This script will filter output from HybPiper based on the output of "hybpiper retrieve_sequences".

As of HybPiper version 2.1.6, "hybpiper retrieve_sequences" only supports filtering based
on project-wide thresholds (i.e. number of total genes recovered). This script will allow
filtering based on individual genes and the mean length or minimum length threshold.

1. Run "hybpiper stats" to generate the "stats.tsv" and "lengths.tsv" files.
2. Run "hybpiper retrieve_sequences "to create a folder of FASTA sequences.
3. Run this command to create new FASTA files based on the per-gene filters provided.
   Also writes the per-gene "deny list" text file.

The FASTA sequences are expected to have the naming scheme of HybPiper:
    geneName.FNA for nucleotide exon files.
    geneName.FAA for amino acid files.
    geneName_supercontig.fasta for supercontig files.
    geneName_intron.fasta for intron-only files.

The geneNames will be taken from either the "hybpiper stats" file (supplied via parameter "--seq_lengths_file")
or a list of gene/sample combinations (supplied via parameter "--denylist"; also produced by running this script).

If you wish to filter intron or supercontig sequences, run this command again with the parameter "--denylist" to
skip the filtering based on lengths.
"""


import os
import sys
import glob
import argparse
import logging
import textwrap
from Bio import SeqIO
from hybpiper import utils
from hybpiper.version import __version__

# Create a custom logger

# Log to Terminal (stderr):
console_handler = logging.StreamHandler(sys.stderr)
console_handler.setLevel(logging.INFO)

# Setup logger:
logger = logging.getLogger(f'hybpiper.{__name__}')

# Add handlers to the logger
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)  # Default level is 'WARNING'

# Get current working directory:
cwd = os.getcwd()


def filter_fastas(deny_dict,
                  sequence_type,
                  sequence_dir,
                  filtered_dir):
    """

    :param dict deny_dict:
    :param str sequence_type:
    :param str sequence_dir: path to directory containing fasta sequence files; default is None
    :param str filtered_dir: directory name for output files, default is None
    """

    # Search within a user-supplied directory for the given fasta sequence files, or the current directory if not:
    if sequence_dir:
        logger.info(f'{"[INFO]:":10} Parsing inpout sequences in directory: {sequence_dir}')
    else:
        sequence_dir = cwd
        logger.info(f'{"[INFO]:":10} Parsing inpout sequences in directory: {cwd}')

    # Set output directory:
    if filtered_dir:
        filtered_dir = filtered_dir
        if not os.path.isdir(filtered_dir):
            os.mkdir(filtered_dir)
    else:
        filtered_dir = '.'

    # Get the correct filename extension:
    seq_suffix = None

    if sequence_type == "supercontig":
        seq_suffix = "_supercontig.fasta"
    elif sequence_type == "intron":
        seq_suffix = "_intron.fasta"
    elif sequence_type == "dna":
        seq_suffix = ".FNA"
    elif sequence_type == "aa":
        seq_suffix = ".FAA"

    # Iterate over the gene fasta files and filter; write filtered output to file:
    fastafiles = list(glob.glob(f'{sequence_dir}/*{seq_suffix}'))

    if len(fastafiles) == 0:
        logger.error(f'{"[ERROR]:":10} There are no sequence files with suffix "{seq_suffix}" in directory:'
                     f' "{sequence_dir}"')
        sys.exit(1)

    for f in fastafiles:
        f_basename = os.path.basename(f)
        gene_name = f_basename.replace(seq_suffix, '')
        genedenylist = set(deny_dict[gene_name])
        new_fn = f'{filtered_dir}/{gene_name}.filtered{seq_suffix}'

        with open(new_fn, 'w') as outfile:
            for seq in SeqIO.parse(f, 'fasta'):
                if seq.id in genedenylist:
                    continue
                else:
                    SeqIO.write(seq, outfile, 'fasta')
    return


def write_denylist(deny_dict,
                   outfile='denylist.txt'):
    """
    Write a text file deny list

    :param dict deny_dict:
    :param str outfile: filename for the "deny list" text file
    """

    with open(outfile, 'w') as outfile_handle:
        for gene in sorted(deny_dict):
            samples = ",".join(deny_dict[gene])
            outfile_handle.write(f'{gene}\t{samples}\n')
    logger.info(f'{"[INFO]:":10} A "deny list" text file has been written to file: {outfile}')

    return


def filter_seqs(gene_lengths_dict,
                minlength,
                minpercent):
    """
    Takes the sample-gene lengths and filters and returns a dictionary by gene of samples to be on the denylist

    :param dict gene_lengths_dict:
    :param int minlength:
    :param float minpercent:
    """

    deny_dict = {}
    total_deny = 0

    for gene in gene_lengths_dict:
        deny_dict[gene] = []
        percent_thresh = gene_lengths_dict[gene]["mean_length"] * minpercent

        for sample_name in gene_lengths_dict[gene]["sample_lengths"]:
            sample_length = gene_lengths_dict[gene]["sample_lengths"][sample_name]

            if sample_length == 0.0:
                continue

            if sample_length < minlength:
                deny_dict[gene].append(sample_name)
                total_deny += 1
                continue

            if sample_length < percent_thresh:
                deny_dict[gene].append(sample_name)
                total_deny += 1

    logger.info(f'{"[INFO]:":10} Filtered out {total_deny} total sequences for {len(deny_dict)} genes based on the'
                f'parameters provided.')

    return deny_dict


def parse_seqlengths(seq_lengths_file):
    """
    Takes the file name for the seqlengths output of hybpiper stats and returns:

    - a list of sample names
    - a dictionary for each gene containing:
        * the name of the gene as the dict key
        * "mean length": integer
        * "sample_lengths": {a dictionary of key:sample_lengths}

    :param str seq_lengths_file: path to the seq_lentghs file output by "hybpiper stats"
    """

    sample_names = []
    gene_lengths_dict = {}

    with open(seq_lengths_file) as seq_lengths_handle:

        gene_names = seq_lengths_handle.readline().rstrip().split("\t")[1:]
        mean_lengths = seq_lengths_handle.readline().rstrip().split("\t")[1:]

        for gene_num in range(len(gene_names)):
            gene_lengths_dict[gene_names[gene_num]] = {"mean_length": float(mean_lengths[gene_num]), "sample_lengths": {}}

        for line in seq_lengths_handle:
            line = line.rstrip().split("\t")
            sample_name = line.pop(0)
            sample_names.append(sample_name)

            for gene_num in range(len(gene_names)):
                gene_lengths_dict[gene_names[gene_num]]["sample_lengths"][sample_name] = float(line[gene_num])

    return sample_names, gene_lengths_dict


def parse_denylist(denylist_fn):
    """
    Parses the text file at denylist_fn and returns a dict with the geneName:[samplelist] pairs

    :param str denylist_fn: path to the text file containing gene/samples combinations to filter out
    """

    if not utils.file_exists_and_not_empty(denylist_fn):
        logger.error(f'{"[ERROR]:":10} File {denylist_fn} does not exist or is empty!')
        sys.exit(1)
    else:
        logger.info(f'{"[INFO]:":10} Parsing "deny list" text file: {denylist_fn}')

    deny_dict = {}
    total_deny = 0

    for line in open(denylist_fn):
        line = line.rstrip().split("\t")
        try:
            samples = line[1].split(",")
        except IndexError:
            samples = []

        total_deny += len(samples)
        deny_dict[line[0]] = samples

    logger.info(f'{"[INFO]:":10} Found {total_deny} total samples at {len(deny_dict)} genes')
    return deny_dict


def standalone():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    group_1 = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('sequence_type',
                        choices=["dna", "aa", "supercontig", "intron"],
                        help='File sequence type for all FASTA files to filter in current directory. For example, the '
                             'amino-acid output of HybPiper would be: aa')
    group_1.add_argument('--denylist',
                         default=None,
                         help='Text file containing gene-sample combinations to omit.\nThe format of the file should '
                              'be one gene per line, a tab, and then a comma-delimited list of samples to disallow: '
                              '\n\n\tgene[tab]sample,sample,sample ')
    group_1.add_argument('--seq_lengths_file',
                         help='Filename for the seq_lengths file (output of the "hybpiper stats" command), with a list '
                              'of genes in the first row, mean target lengths in the second row, and sample recovery '
                              'in other rows.')
    parser.add_argument('--denylist_filename',
                        default='denylist.txt',
                        type=str,
                        help='File name for the "deny list" text file (if written). Default is <denylist.txt>')
    parser.add_argument('--length_filter',
                        default=0,
                        type=int,
                        help='Minimum length to allow a sequence in nucleotides for DNA or amino acids for protein '
                             'sequences')
    parser.add_argument('--percent_filter',
                        default=0,
                        type=float,
                        help='Minimum fraction (between 0 and 1) of the mean target length to allow a sequence for a '
                             'gene. Lengths taken from HybPiper stats file.')
    parser.add_argument('--sequence_dir',
                        default=None,
                        help='Specify directory containing sequences output by the "hybpiper retrieve_sequences" '
                             'command. Default is to search in the current working directory')
    parser.add_argument('--filtered_dir',
                        default=None,
                        help='Specify directory for output filtered FASTA files. Default is to write to the current '
                             'working directory')

    args = parser.parse_args()
    main(args)


def main(args):
    """
    Entry point for the assemble.py module.

    :param argparse.Namespace args:
    """

    logger.info(f'{"[INFO]:":10} HybPiper version {__version__} was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.info(f'{fill}\n')

    # Check some args:
    if args.denylist:

        if args.sequence_type not in ["supercontig", "intron"]:
            logger.error(f'{"[ERROR]:":10} You have provided a file via the flag "--denylist". Please use '
                         f'"--sequence_type intron" or "--sequence_type supercontig"')
            sys.exit(1)
        else:
            deny_dict = parse_denylist(args.denylist)
    else:
        sample_names, gene_lengths_dict = parse_seqlengths(args.seq_lengths_file)

        deny_dict = filter_seqs(gene_lengths_dict,
                                args.length_filter,
                                args.percent_filter)

        write_denylist(deny_dict,
                       outfile=args.denylist_filename)

    # Filter the fasta sequences based on the deny_dict:
    filter_fastas(deny_dict,
                  args.sequence_type,
                  args.sequence_dir,
                  args.filtered_dir)


########################################################################################################################
########################################################################################################################

if __name__ == '__main__':
    standalone()

########################################################################################################################
########################################################################################################################
