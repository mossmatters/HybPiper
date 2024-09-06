#!/usr/bin/env python

#########################
# HybPiper Stats Script #
#########################

"""
Writes a report file called "seq_lengths.tsv" (default filename, user can change this). The first line contains gene
names. The second line contains the length of the reference sequences (targets). If there are multiple targets per gene,
the mean length is reported. All other rows contain one sample per line.

Parses the "seq_lengths.tsv" and gathers additional statistics about the HybPiper run. Writes a report file called
"hybpiper_stats.tsv" (default filename, user can change this).

For an explanation of columns, see github.com/mossmatters/HybPiper/wiki
"""

import argparse
import os
import sys
import subprocess
import re
from Bio import SeqIO
from collections import defaultdict
import logging
import textwrap

from hybpiper import utils
from hybpiper.version import __version__


# Create a custom logger

# Log to Terminal (stderr):
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)

# Setup logger:
logger = logging.getLogger(f'hybpiper.{__name__}')

# Add handlers to the logger
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)  # Default level is 'WARNING'


def get_seq_lengths(targetfile,
                    targetfile_sequence_type,
                    list_of_sample_names,
                    sequence_type_to_calculate_stats_for,
                    seq_lengths_filename,
                    compressed_sample_dict,
                    sampledir_parent):
    """
    Recover the sequence length of each target file gene (calculated as mean length if a representative sequence from
    more than one taxon is provided for a given gene). Calculate the percentage length recovery for each gene,
    for each sample. If a protein target file was used, convert target gene lengths to the number of nucleotides
    (i.e. amino-acids x 3).

    :param str targetfile: path to the targetfile
    :param str targetfile_sequence_type: sequence type in the target file ('DNA' or 'protein')
    :param list list_of_sample_names: a list of sample names
    :param str sequence_type_to_calculate_stats_for: gene (in nucleotides) or supercontig (in nucleotides)
    :param str seq_lengths_filename: optional filename for seq_lengths file. Default is seq_lengths.tsv
    :param dict compressed_sample_dict: dictionary containing *.tar.gz contents for compressed sample directories
    :param path sampledir_parent: path to the parent directory containing the sample directories
    :return str seq_lengths_report_filename: path to the sequence length report file written by this function
    """

    lines_for_report = []  # lines to write to file

    # Set variable 'filetype', used to reconstruct path to each sequence:
    filetype = None
    if sequence_type_to_calculate_stats_for.upper() == 'GENE':
        filetype = 'FNA'
    elif sequence_type_to_calculate_stats_for.upper() == "SUPERCONTIG":
        filetype = 'supercontig'
    assert filetype

    # Get the names and lengths for each sequence in the target file:
    gene_names = []
    reference_lengths = defaultdict(list)
    for prot in SeqIO.parse(targetfile, "fasta"):
        protname = prot.id.split("-")[-1]
        gene_names.append(protname)
        if targetfile_sequence_type.upper() == 'PROTEIN':
            reference_lengths[protname].append(len(prot.seq) * 3)  # covert from amino-acids to nucleotides
        elif targetfile_sequence_type.upper() == 'DNA':
            reference_lengths[protname].append(len(prot.seq))

    unique_names = list(set(gene_names))
    avg_ref_lengths = [(sum(reference_lengths[gene])/len(reference_lengths[gene])) for gene in unique_names]

    # Capture the unique gene names and average length for each gene to write to a report:
    unique_names_to_write = '\t'.join(unique_names)
    avg_ref_lengths_to_write = '\t'.join([str(x) for x in avg_ref_lengths])
    lines_for_report.append(f'Species\t{unique_names_to_write}')
    lines_for_report.append(f'MeanLength\t{avg_ref_lengths_to_write}')

    # Get seq lengths for sample gene sequences (FNA or supercontigs):
    sample_name_to_total_bases_dict = defaultdict(int)

    for sample_name in list_of_sample_names:  # iterate over sample names

        name_lengths = []  # lengths of sequences in nucleotides
        for gene in range(len(unique_names)):  # iterate over genes

            # Reconstruct path to the sequence with sample directory as root:
            if filetype == 'supercontig':
                seq_file = os.path.join(sample_name, unique_names[gene], sample_name, 'sequences', 'intron',
                                        f'{unique_names[gene]}_supercontig.fasta')
            else:
                seq_file = os.path.join(sample_name, unique_names[gene], sample_name, "sequences", filetype,
                                        f'{unique_names[gene]}.{filetype}')

            # Get full path to seq_file:
            seq_file_full_path = f'{sampledir_parent}/{seq_file}'

            # Get sequence details depending on compressed vs non-compressed sample directory:
            seq_length = None

            if sample_name in compressed_sample_dict:
                seq_file_exists = True if seq_file in compressed_sample_dict[sample_name] else False
                seq_file_size = compressed_sample_dict[sample_name][seq_file] \
                    if seq_file in compressed_sample_dict[sample_name] else 0

                if seq_file_exists:
                    if seq_file_size == 0:
                        logger.warning(f'{"[WARNING]:":10} File {seq_file_full_path} exists, but is empty! A length of '
                                       f'zero will be recorded for this sequence.\n')
                        name_lengths.append("0")
                    else:
                        seqrecord = utils.get_compressed_seqrecord(sample_name,
                                                                   sampledir_parent,
                                                                   seq_file)
                        seq_length = len(seqrecord.seq.replace('N', ''))
                        if seq_length > 1.5 * avg_ref_lengths[gene] and filetype != 'supercontig':
                            logger.warning(f'{"[WARNING]:":10} Sequence length for {sample_name} is more than 50% '
                                           f'longer than {unique_names[gene]} reference!\n')

                        name_lengths.append(str(seq_length))
                        sample_name_to_total_bases_dict[sample_name] += seq_length
                else:
                    name_lengths.append("0")

            else:  # i.e. uncompressed sample folder

                seq_file_exists = True if os.path.isfile(seq_file_full_path) else False
                seq_file_size = os.path.getsize(seq_file_full_path) if os.path.isfile(seq_file_full_path) else 0

                if seq_file_exists:
                    if seq_file_size == 0:
                        logger.warning(f'{"[WARNING]:":10} File {seq_file_full_path} exists, but is empty! A length of '
                                       f'zero will be recorded for this sequence.\n')
                        name_lengths.append("0")
                    else:
                        seq_length = len(SeqIO.read(seq_file_full_path, 'fasta').seq.replace('N', ''))
                        if seq_length > 1.5 * avg_ref_lengths[gene] and filetype != 'supercontig':
                            logger.warning(f'{"[WARNING]:":10} Sequence length for {sample_name} is more than 50% '
                                           f'longer than {unique_names[gene]} reference!\n')

                        name_lengths.append(str(seq_length))
                        sample_name_to_total_bases_dict[sample_name] += seq_length
                else:
                    name_lengths.append("0")

        lengths_to_write = '\t'.join(name_lengths)
        lines_for_report.append(f'{sample_name}\t{lengths_to_write}')

    # Write report file "seq_lengths.tsv"
    seq_lengths_report_filename = f'{seq_lengths_filename}.tsv'
    with open(seq_lengths_report_filename, 'w') as seq_lengths_handle:
        for item in lines_for_report:
            seq_lengths_handle.write(f'{item}\n')

    logger.info(f'{"[INFO]:":10} A sequence length table has been written to file: {seq_lengths_filename}.tsv')

    return seq_lengths_report_filename, sample_name_to_total_bases_dict


def file_len(fname):
    """
    Function to recover the number of lines in a test file. Runs the command-line builtin 'wc -l' via subprocess.

    :param str fname: path to a file (spades_genelist.txt/genes_with_seqs.txt/exonerate_genelist.txt)
    :return int: number of lines in the file as reported by the command-line builtin 'wc -l'
    """

    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])


def enrich_efficiency_blastx(sample_name,
                             sampledir_parent,
                             blastxfilename,
                             compressed_sample_bool=False,
                             blastxfile_unpaired_exists=False,
                             total_input_reads_paired_exists=False,
                             total_input_reads_single_exists=False,
                             total_input_reads_unpaired_exists=False):
    """
    Parse BLASTX results to calculate enrichment efficiency

    :param str sample_name: name of a given sample
    :param path sampledir_parent:
    :param str blastxfilename: path to the *.tsv BLASTx output filename for a given sample
    :param bool compressed_sample_bool: True if sample directory is a compressed tarball
    :param bool blastxfile_unpaired_exists: True if an unpaired BLASTX exists for this sample
    :param bool total_input_reads_paired_exists:
    :param bool total_input_reads_single_exists:
    :param bool total_input_reads_unpaired_exists:
    :return str, str, str: values for input reads, mapped reads, and percent mapped reads
    """

    # Set expected file paths with sample folder as root:
    unpaired_blastx_file = blastxfilename.replace(".blastx", "_unpaired.blastx")
    total_input_reads_paired = f'{sample_name}/total_input_reads_paired.txt'
    total_input_reads_single = f'{sample_name}/total_input_reads_single.txt'
    total_input_reads_unpaired = f'{sample_name}/total_input_reads_unpaired.txt'

    # Process stats from the 'main' BLAST file:
    if compressed_sample_bool:
        blastx_lines = utils.get_compressed_file_lines(sample_name,
                                                       sampledir_parent,
                                                       blastxfilename)
    else:
        blastx_handle = open(f'{sampledir_parent}/{blastxfilename}', 'r')
        blastx_lines = blastx_handle.readlines()

    reads_with_hits = [x.split('\t')[0] for x in blastx_lines if x]

    if not compressed_sample_bool:
        blastx_handle.close()

    # If an unpaired BLASTX exists for this sample, process stats:
    if blastxfile_unpaired_exists:
        if compressed_sample_bool:
            blastx_unpaired_lines = utils.get_compressed_file_lines(sample_name,
                                                                    sampledir_parent,
                                                                    unpaired_blastx_file)
        else:
            blastx_unpaired_handle = open(f'{sampledir_parent}/{unpaired_blastx_file}', 'r')
            blastx_unpaired_lines = blastx_unpaired_handle.readlines()

        reads_with_hits += [x.split('\t')[0] for x in blastx_unpaired_lines if x]

        if not compressed_sample_bool:
            blastx_unpaired_handle.close()

    mapped_reads = len(set(reads_with_hits))

    # Recover total numbers of reads in input files:
    if total_input_reads_paired_exists:

        if compressed_sample_bool:
            total_input_reads_paired_lines = utils.get_compressed_file_lines(sample_name,
                                                                             sampledir_parent,
                                                                             total_input_reads_paired)
            total_input_reads = int(total_input_reads_paired_lines[0].rstrip())
        else:
            with open(f'{sampledir_parent}/{total_input_reads_paired}', 'r') as paired_number:
                total_input_reads = int(paired_number.read().rstrip())

    elif total_input_reads_single_exists:

        if compressed_sample_bool:
            total_input_reads_single_lines = utils.get_compressed_file_lines(sample_name,
                                                                             sampledir_parent,
                                                                             total_input_reads_single)
            total_input_reads = int(total_input_reads_single_lines[0].rstrip())
        else:
            with open(f'{sampledir_parent}/{total_input_reads_single}', 'r') as single_number:
                total_input_reads = int(single_number.read().rstrip())
    else:
        fill = utils.fill_forward_slash(
            f'{"[WARNING]:":10} No file containing total input paired or single-end read count found for '
            f'sample {sample_name}!. Please check the log file for this sample. If "hybpiper assemble" can not be '
            f'run successfully for this sample, please remove it from your namelist text file before running '
            f'"hybpiper stats".',
            width=90, subsequent_indent=' ' * 11, break_long_words=False,
            break_on_forward_slash=True)

        logger.warning(f'{fill}')
        sys.exit()

    if blastxfile_unpaired_exists and total_input_reads_unpaired_exists:

        if compressed_sample_bool:
            total_input_reads_unpaired_lines = utils.get_compressed_file_lines(sample_name,
                                                                               sampledir_parent,
                                                                               total_input_reads_unpaired)

            total_input_reads = total_input_reads + int(total_input_reads_unpaired_lines[0].rstrip())
        else:
            with open(f'{sampledir_parent}/{total_input_reads_unpaired}', 'r') as unpaired_number:
                total_input_reads = total_input_reads + int(unpaired_number.read().rstrip())

    try:
        pct_mapped = 100 * mapped_reads / total_input_reads
    except ZeroDivisionError:
        pct_mapped = 0.0

    return str(total_input_reads), str(mapped_reads), "{0:.1f}".format(pct_mapped)


def enrich_efficiency_bwa(sample_name,
                          sampledir_parent,
                          compressed_sample_dict,
                          compressed_sample_bool=False,
                          bam_file_unpaired_exists=False):
    """
    Run and parse samtools flagstat output, return number of reads and number on target. Calculate percentage of
    reads mapped.

    :param str sample_name: sample name
    :param path sampledir_parent: sampledir_parent name
    :param dict compressed_sample_dict: dictionary containing *.tar.gz contents for compressed sample directories
    :param bool compressed_sample_bool: True if sample directory is a compressed tarball
    :param bool bam_file_unpaired_exists: True if an unpaired bamfile exists for this sample
    :return str, str, str: values for input reads, mapped reads, and percent mapped reads:
    """

    # Set expected file paths with sample folder as root:
    bam_file = f'{sample_name}/{sample_name}.bam'
    bam_file_unpaired = f'{sample_name}/{sample_name}_unpaired.bam'
    bam_flagstats_tsv_file = f'{sample_name}/{sample_name}_bam_flagstat.tsv'
    unpaired_bam_flagstats_tsv_file = f'{sample_name}/{sample_name}_unpaired_bam_flagstat.tsv'

    ####################################################################################################################
    # Check if the bam flagstat file(s) exist (i.e., samples were run with HybPiper >= v2.3.0). If not, warn the user
    # and exit:
    ####################################################################################################################
    if compressed_sample_bool:

        bam_flagstats_tsv_file_exists = \
            True if bam_flagstats_tsv_file in compressed_sample_dict[sample_name] else False

        unpaired_bam_flagstats_tsv_file_exists = \
            True if unpaired_bam_flagstats_tsv_file in compressed_sample_dict[sample_name] else False

    else:
        bam_flagstats_tsv_file_exists = \
            True if utils.file_exists_and_not_empty(f'{sampledir_parent}/{bam_flagstats_tsv_file}') else False

        unpaired_bam_flagstats_tsv_file_exists = \
            True if utils.file_exists_and_not_empty(f'{sampledir_parent}/{unpaired_bam_flagstats_tsv_file}') else False

    ####################################################################################################################
    # Process stats from the 'main' bam file:
    ####################################################################################################################

    # Initialise count at zero:
    num_reads = 0
    mapped_reads = 0

    if compressed_sample_bool:
        if bam_flagstats_tsv_file_exists:

            bam_flagstat_lines = utils.get_compressed_file_lines(sample_name,
                                                                 sampledir_parent,
                                                                 bam_flagstats_tsv_file)
        else:
            bam_flagstat_lines = utils.get_bamtools_flagstat_lines_from_compressed(sample_name,
                                                                                   sampledir_parent,
                                                                                   bam_file)
    else:
        if bam_flagstats_tsv_file_exists:

            bam_flagstat_handle = open(f'{sampledir_parent}/{bam_flagstats_tsv_file}', 'r')
            bam_flagstat_lines = list(bam_flagstat_handle.readlines())
            bam_flagstat_handle.close()

        else:
            bam_flagstat_lines = utils.get_bamtools_flagstat_lines_from_uncompressed(sample_name,
                                                                                     sampledir_parent,
                                                                                     bam_file)

    for line in bam_flagstat_lines:
        if re.search('primary$', line):
            num_reads = float(line.split('\t')[0])
        if re.search(r'\bprimary mapped$', line):
            mapped_reads = float(line.split('\t')[0])

    ####################################################################################################################
    # If an unpaired bam exists for this sample, process stats:
    ####################################################################################################################
    if bam_file_unpaired_exists:

        if compressed_sample_bool:

            if unpaired_bam_flagstats_tsv_file_exists:
                bam_flagstat_unpaired_lines = utils.get_compressed_file_lines(sample_name,
                                                                              sampledir_parent,
                                                                              unpaired_bam_flagstats_tsv_file)
            else:
                bam_flagstat_unpaired_lines = \
                    utils.get_bamtools_flagstat_lines_from_compressed(sample_name,
                                                                      sampledir_parent,
                                                                      bam_file_unpaired)
        else:

            if unpaired_bam_flagstats_tsv_file_exists:

                bam_flagstat_unpaired_handle = open(f'{sampledir_parent}/{unpaired_bam_flagstats_tsv_file}', 'r')
                bam_flagstat_unpaired_lines = list(bam_flagstat_unpaired_handle.readlines())
                bam_flagstat_unpaired_handle.close()

            else:
                bam_flagstat_unpaired_lines = \
                    utils.get_bamtools_flagstat_lines_from_uncompressed(sample_name,
                                                                        sampledir_parent,
                                                                        bam_file_unpaired)

        for line in bam_flagstat_unpaired_lines:
            if re.search('primary$', line):
                num_reads += float(line.split('\t')[0])
            if re.search(r'\bprimary mapped\b', line):
                mapped_reads += float(line.split('\t')[0])

    try:
        pct_mapped = 100 * mapped_reads / num_reads
    except ZeroDivisionError:
        pct_mapped = 0.0

    return str(int(num_reads)), str(int(mapped_reads)), "{0:.1f}".format(pct_mapped)


def recovery_efficiency(sample_name,
                        sampledir_parent,
                        compressed_sample_bool,
                        compressed_sample_dict):
    """
    Reports the number of genes with mapping hits, contigs, and exon sequences

    :param str sample_name: sample name
    :param str sampledir_parent: sampledir_parent name
    :param bool compressed_sample_bool: True if sample directory is a compressed tarball
    :param dict compressed_sample_dict: dictionary containing *.tar.gz contents for compressed sample directories
    :return list: a list containing the number of genes with contigs, exonerate hits, and assembled sequences
    """

    txt_files = [
                 "spades_genelist.txt",
                 "exonerate_genelist.txt",
                 "genes_with_seqs.txt"
                 ]

    my_stats = []
    for txt in txt_files:
        sample_file_path = f'{sample_name}/{txt}'  # sample folder as root

        if compressed_sample_bool:
            if sample_file_path in compressed_sample_dict[sample_name]:

                sample_file_path_lines = utils.get_compressed_file_lines(sample_name,
                                                                         sampledir_parent,
                                                                         sample_file_path)
                number_of_lines = len(sample_file_path_lines)
                my_stats.append(number_of_lines)
            else:
                my_stats.append(0)
        else:
            if os.path.isfile(f'{sampledir_parent}/{sample_file_path}'):
                my_stats.append(file_len(f'{sampledir_parent}/{sample_file_path}'))
            else:
                my_stats.append(0)

    return [str(a) for a in my_stats]


def seq_length_calc(seq_lengths_fn):
    """
    From the output file produced by get_seq_lengths(), calculate the number of genes with seqs, and at least a
    percentage of the reference length.

    :param path seq_lengths_fn: path to the sequence length file produced by get_seq_lengths()
    :return dict seq_length_dict: dictionary of sample:list_of_genes_above_length_threshold
    """

    seq_length_dict = {}

    with open(seq_lengths_fn) as seq_len:
        gene_names = seq_len.readline()  # skip the first line
        target_lengths = seq_len.readline().split()[1:]

        for line in seq_len:
            line = line.split()
            name = line.pop(0)
            is_25pct = 0
            is_50pct = 0
            is_75pct = 0
            is_150pct = 0

            for gene in range(len(line)):
                gene_length = float(line[gene])
                target_length = float(target_lengths[gene])
                if gene_length > target_length * 0.25:
                    is_25pct += 1
                if gene_length > target_length * 0.50:
                    is_50pct += 1
                if gene_length > target_length * 0.75:
                    is_75pct += 1
                if gene_length > target_length * 1.5:
                    is_150pct += 1

            seq_length_dict[name] = [str(is_25pct), str(is_50pct), str(is_75pct), str(is_150pct)]

    return seq_length_dict


def standalone():
    """
    Used when this module is run as a stand-alone script. Parses command line arguments and runs function main().:
    """

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    group_1 = parser.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna', dest='targetfile_dna',
                         help='FASTA file containing DNA target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa', dest='targetfile_aa',
                         help='FASTA file containing amino-acid target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    parser.add_argument("sequence_type",
                        help="Sequence type (gene or supercontig) to recover stats for",
                        choices=["gene", "GENE", "supercontig", "SUPERCONTIG"])
    parser.add_argument("namelist",
                        help="text file with names of HybPiper output directories, one per line")
    parser.add_argument("--seq_lengths_filename",
                        help="File name for the sequence lengths *.tsv file. Default is <seq_lengths.tsv>.",
                        default='seq_lengths')
    parser.add_argument("--stats_filename",
                        help="File name for the stats *.tsv file. Default is <hybpiper_stats.tsv>",
                        default='hybpiper_stats')

    parser.set_defaults(targetfile_dna=False, targetfile_aa=False)

    args = parser.parse_args()
    main(args)


def main(args):
    """
    Entry point for the hybpiper_main.py module.

    :param argparse.Namespace args:
    """

    logger.info(f'{"[INFO]:":10} HybPiper version {__version__} was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.info(f'{fill}\n')

    logger.info(f'{"[INFO]:":10} Recovering statistics for the HybPiper run(s)...')

    ####################################################################################################################
    # Set target file name and type:
    ####################################################################################################################
    targetfile = args.targetfile_dna if args.targetfile_dna else args.targetfile_aa
    targetfile_type = 'DNA' if args.targetfile_dna else 'protein'

    ####################################################################################################################
    # Check for presence of required input files:
    ####################################################################################################################
    logger.info(f'{"[INFO]:":10} The following file of sample names was provided: "{args.namelist}".')
    if not utils.file_exists_and_not_empty(args.namelist):
        sys.exit(f'{"[ERROR]:":10} File {args.namelist} is missing or empty, exiting!')

    logger.info(f'{"[INFO]:":10} The following target file was provided: "{targetfile}".')
    if not utils.file_exists_and_not_empty(targetfile):
        sys.exit(f'{"[ERROR]:":10} File {targetfile} is missing or empty, exiting!')

    ####################################################################################################################
    # Search within a user-supplied directory for the given sample directories, or the current directory if not:
    ####################################################################################################################
    if args.hybpiper_dir:
        sampledir_parent = os.path.abspath(args.hybpiper_dir)
    else:
        sampledir_parent = os.getcwd()

    if not os.path.isdir(sampledir_parent):
        sys.exit(f'{"[ERROR]:":10} Folder {sampledir_parent} not found, exiting!')

    ####################################################################################################################
    # Parse namelist and check for the presence of the corresponding sample directories or *.tar.gz files:
    ####################################################################################################################
    list_of_sample_names = []

    with open(args.namelist, 'r') as namelist_handle:
        for line in namelist_handle.readlines():
            sample_name = line.rstrip()
            if sample_name:
                list_of_sample_names.append(sample_name)
                if re.search('/', sample_name):
                    sys.exit(f'{"[ERROR]:":10} A sample name must not contain '
                             f'forward slashes. The file {args.namelist} contains: {sample_name}')

    samples_found = []
    compressed_sample_dict = dict()

    for sample_name in list_of_sample_names:
        compressed_sample = f'{sampledir_parent}/{sample_name}.tar.gz'
        uncompressed_sample = f'{sampledir_parent}/{sample_name}'

        if os.path.isfile(compressed_sample):
            compressed_sample_dict[sample_name] = utils.parse_compressed_sample(compressed_sample)
            samples_found.append(sample_name)

        if os.path.isdir(uncompressed_sample):
            if sample_name in samples_found:
                sys.exit(f'{"[ERROR]:":10} Both a compressed and an un-compressed sample folder have been found for '
                         f'sample {sample_name} in directory {sampledir_parent}. Please remove one!')
            else:
                samples_found.append(sample_name)

    samples_missing = sorted(list(set(list_of_sample_names) - set(samples_found)))

    if samples_missing:

        fill = utils.fill_forward_slash(f'{"[WARNING]:":10} File {args.namelist} contains samples not found in '
                                        f'directory "{sampledir_parent}". The missing samples are:',
                                        width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                        break_on_forward_slash=True)

        logger.warning(f'{fill}\n')

        for name in samples_missing:
            logger.warning(f'{" " * 10} {name}')
        logger.warning('')

    list_of_sample_names = [x for x in list_of_sample_names if x not in samples_missing]

    ####################################################################################################################
    # Get sequence lengths for recovered genes, and write them to file, along with total bases recovered:
    ####################################################################################################################
    (seq_lengths_file_path,
     sample_name_to_total_bases_dict) = get_seq_lengths(targetfile,
                                                        targetfile_type,
                                                        list_of_sample_names,
                                                        args.sequence_type,
                                                        args.seq_lengths_filename,
                                                        compressed_sample_dict,
                                                        sampledir_parent)

    seq_length_dict = seq_length_calc(seq_lengths_file_path)

    lines_for_stats_report = []

    categories = [
                  "Name",
                  "NumReads",
                  "ReadsMapped",
                  "PctOnTarget",
                  "GenesMapped",
                  "GenesWithContigs",
                  "GenesWithSeqs",
                  "GenesAt25pct",
                  "GenesAt50pct",
                  "GenesAt75pct",
                  "GenesAt150pct",
                  "ParalogWarningsLong",
                  "ParalogWarningsDepth",
                  "GenesWithoutStitchedContigs",
                  "GenesWithStitchedContigs",
                  "GenesWithStitchedContigsSkipped",
                  "GenesWithChimeraWarning",
                  "TotalBasesRecovered"
                  ]

    categories_for_printing = '\t'.join(categories)
    lines_for_stats_report.append(categories_for_printing)

    ####################################################################################################################
    # Iterate over sample names and populate stats_dict:
    ####################################################################################################################
    stats_dict = {}
    for sample_name in list_of_sample_names:

        # Check is the sample directory is a compressed tarball:
        compressed_sample_bool = True if sample_name in compressed_sample_dict else False

        # Initialise stats list for the sample:
        stats_dict[sample_name] = []

        ################################################################################################################
        # Enrichment efficiency:
        ################################################################################################################

        # Set expected paths with sample folder as root:
        bamfile = f'{sample_name}/{sample_name}.bam'
        bamfile_unpaired = f'{sample_name}/{sample_name}_unpaired.bam'
        blastxfile = f'{sample_name}/{sample_name}.blastx'
        blastxfile_unpaired = f'{sample_name}/{sample_name}_unpaired.blastx'
        total_input_reads_paired = f'{sample_name}/total_input_reads_paired.txt'
        total_input_reads_single = f'{sample_name}/total_input_reads_single.txt'
        total_input_reads_unpaired = f'{sample_name}/total_input_reads_unpaired.txt'

        if compressed_sample_bool:

            # Check if files are present in the compressed folder:
            bam_file_exists = True if bamfile in compressed_sample_dict[sample_name] else False
            bam_file_unpaired_exists = True if bamfile_unpaired in compressed_sample_dict[sample_name] else False
            blastxfile_exists = True if blastxfile in compressed_sample_dict[sample_name] else False
            blastxfile_unpaired_exists = True if blastxfile_unpaired in compressed_sample_dict[sample_name] else False

            total_input_reads_paired_exists = \
                True if total_input_reads_paired in compressed_sample_dict[sample_name] else False
            total_input_reads_single_exists = \
                True if total_input_reads_single in compressed_sample_dict[sample_name] else False
            total_input_reads_unpaired_exists = \
                True if total_input_reads_unpaired in compressed_sample_dict[sample_name] else False

            if not bam_file_exists and not blastxfile_exists:
                fill = utils.fill_forward_slash(f'{"[WARNING]:":10} No *.bam or *.blastx file found for '
                                                f'{sample_name}. No statistics will be recovered. Please check the log '
                                                f'file for this sample!',
                                                width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                                break_on_forward_slash=True)
                logger.warning(f'{fill}')
                continue

            if bam_file_exists:
                stats_dict[sample_name] += enrich_efficiency_bwa(
                    sample_name,
                    sampledir_parent,
                    compressed_sample_dict,
                    compressed_sample_bool=compressed_sample_bool,
                    bam_file_unpaired_exists=bam_file_unpaired_exists
                )

            else:
                stats_dict[sample_name] += enrich_efficiency_blastx(
                    sample_name,
                    sampledir_parent,
                    blastxfile,
                    compressed_sample_bool=compressed_sample_bool,
                    blastxfile_unpaired_exists=blastxfile_unpaired_exists,
                    total_input_reads_paired_exists=total_input_reads_paired_exists,
                    total_input_reads_single_exists=total_input_reads_single_exists,
                    total_input_reads_unpaired_exists=total_input_reads_unpaired_exists
                )

        else:  # i.e. uncompressed sample folder

            # Check if files are present in the un-compressed folder:
            bam_file_exists = True if os.path.isfile(f'{sampledir_parent}/{bamfile}') else False
            bam_file_unpaired_exists = True if os.path.isfile(f'{sampledir_parent}/{bamfile_unpaired}') else False
            blastxfile_exists = True if os.path.isfile(f'{sampledir_parent}/{blastxfile}') else False
            blastxfile_unpaired_exists = True if os.path.isfile(f'{sampledir_parent}/{blastxfile_unpaired}') else False

            total_input_reads_paired_exists = \
                True if os.path.isfile(f'{sampledir_parent}/{total_input_reads_paired}') else False
            total_input_reads_single_exists = \
                True if os.path.isfile(f'{sampledir_parent}/{total_input_reads_single}') else False
            total_input_reads_unpaired_exists = \
                True if os.path.isfile(f'{sampledir_parent}/{total_input_reads_unpaired}') else False

            if not bam_file_exists and not blastxfile_exists:
                fill = utils.fill_forward_slash(f'{"[WARNING]:":10} No *.bam or *.blastx file found for '
                                                f'{sample_name}. No statistics will be recovered. Please check the log '
                                                f'file for this sample!',
                                                width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                                break_on_forward_slash=True)
                logger.warning(f'{fill}')
                continue

            if bam_file_exists:
                stats_dict[sample_name] += enrich_efficiency_bwa(sample_name,
                                                                 sampledir_parent,
                                                                 compressed_sample_dict,
                                                                 bam_file_unpaired_exists=bam_file_unpaired_exists)

            else:
                stats_dict[sample_name] += enrich_efficiency_blastx(
                    sample_name,
                    sampledir_parent,
                    blastxfile,
                    blastxfile_unpaired_exists=blastxfile_unpaired_exists,
                    total_input_reads_paired_exists=total_input_reads_paired_exists,
                    total_input_reads_single_exists=total_input_reads_single_exists,
                    total_input_reads_unpaired_exists=total_input_reads_unpaired_exists
                )

        ################################################################################################################
        # Recovery efficiency:
        ################################################################################################################
        stats_dict[sample_name] += recovery_efficiency(sample_name,
                                                       sampledir_parent,
                                                       compressed_sample_bool,
                                                       compressed_sample_dict)

        ################################################################################################################
        # Sequence length:
        ################################################################################################################
        stats_dict[sample_name] += seq_length_dict[sample_name]

        ################################################################################################################
        # Paralogs - long:
        ################################################################################################################
        long_paralog_warnings_file = f'{sample_name}/{sample_name}_genes_with_long_paralog_warnings.txt'

        if compressed_sample_bool:
            if long_paralog_warnings_file in compressed_sample_dict[sample_name]:
                long_paralog_warnings_file_lines = utils.get_compressed_file_lines(sample_name,
                                                                                   sampledir_parent,
                                                                                   long_paralog_warnings_file)
                paralog_warns = len(long_paralog_warnings_file_lines)
                stats_dict[sample_name].append(str(paralog_warns))
            else:
                stats_dict[sample_name].append("0")
        else:
            full_path = f'{sampledir_parent}/{long_paralog_warnings_file}'
            if os.path.isfile(full_path):
                paralog_warns = file_len(full_path)
                stats_dict[sample_name].append(str(paralog_warns))
            else:
                stats_dict[sample_name].append("0")

        ################################################################################################################
        # Paralogs - by contig depth across query protein:
        ################################################################################################################
        depth_paralog_warnings_file = f'{sample_name}/{sample_name}_genes_with_paralog_warnings_by_contig_depth.csv'

        if compressed_sample_bool:
            num_genes_paralog_warning_by_depth = 0

            if depth_paralog_warnings_file in compressed_sample_dict[sample_name]:
                depth_paralog_warnings_file_lines = utils.get_compressed_file_lines(sample_name,
                                                                                    sampledir_parent,
                                                                                    depth_paralog_warnings_file)
                for gene_stats in depth_paralog_warnings_file_lines:
                    stat = gene_stats.split(',')[3].strip()
                    if stat == 'True':
                        num_genes_paralog_warning_by_depth += 1

            stats_dict[sample_name].append(str(num_genes_paralog_warning_by_depth))

        else:
            num_genes_paralog_warning_by_depth = 0
            full_path = f'{sampledir_parent}/{depth_paralog_warnings_file}'
            if os.path.isfile(full_path):
                with open(full_path) as paralogs_by_depth:
                    lines = paralogs_by_depth.readlines()
                    for gene_stats in lines:
                        stat = gene_stats.split(',')[3].strip()
                        if stat == 'True':
                            num_genes_paralog_warning_by_depth += 1

            stats_dict[sample_name].append(str(num_genes_paralog_warning_by_depth))

        ################################################################################################################
        # Stitched contig information:
        ################################################################################################################
        stitched_contig_produced = 0
        no_stitched_contig = 0
        stitched_contig_skipped = 0

        genes_with_stitched_contig_file = f'{sample_name}/{sample_name}_genes_with_stitched_contig.csv'

        if compressed_sample_bool:

            if genes_with_stitched_contig_file in compressed_sample_dict[sample_name]:
                lines = utils.get_compressed_file_lines(sample_name,
                                                        sampledir_parent,
                                                        genes_with_stitched_contig_file)
                for gene_stats in lines:
                    stat = gene_stats.split(',')[2]
                    if re.search('single Exonerate hit', stat):
                        no_stitched_contig += 1
                    elif re.search('Stitched contig produced', stat):
                        stitched_contig_produced += 1
                    elif re.search('Stitched contig step skipped', stat):
                        stitched_contig_skipped += 1

        else:
            full_path = f'{sampledir_parent}/{genes_with_stitched_contig_file}'
            if os.path.isfile(full_path):
                with open(full_path) as stitched_contig_stats:
                    lines = stitched_contig_stats.readlines()
                    for gene_stats in lines:
                        stat = gene_stats.split(',')[2]
                        if re.search('single Exonerate hit', stat):
                            no_stitched_contig += 1
                        elif re.search('Stitched contig produced', stat):
                            stitched_contig_produced += 1
                        elif re.search('Stitched contig step skipped', stat):
                            stitched_contig_skipped += 1

        stitched_contigs_produced_total = stitched_contig_produced
        stats_dict[sample_name].append(str(no_stitched_contig))
        stats_dict[sample_name].append(str(stitched_contigs_produced_total))
        stats_dict[sample_name].append(str(stitched_contig_skipped))

        ################################################################################################################
        # Chimeric stitched contigs:
        ################################################################################################################
        chimeric_stitched_contigs = 0
        genes_derived_from_putative_chimeric_stitched_contig_file = \
            f'{sample_name}/{sample_name}_genes_derived_from_putative_chimeric_stitched_contig.csv'

        if compressed_sample_bool:

            if genes_derived_from_putative_chimeric_stitched_contig_file in compressed_sample_dict[sample_name]:
                lines = utils.get_compressed_file_lines(sample_name,
                                                        sampledir_parent,
                                                        genes_derived_from_putative_chimeric_stitched_contig_file)
                for gene_stats in lines:
                    stat = gene_stats.split(',')[2]
                    if re.search(' Chimera WARNING for stitched contig.', stat):
                        chimeric_stitched_contigs += 1

        else:
            full_path = f'{sampledir_parent}/{genes_derived_from_putative_chimeric_stitched_contig_file}'
            if os.path.isfile(full_path):
                with open(full_path) as chimeric_stitched_contig_stats:

                    lines = chimeric_stitched_contig_stats.readlines()
                    for gene_stats in lines:
                        stat = gene_stats.split(',')[2]
                        if re.search(' Chimera WARNING for stitched contig.', stat):
                            chimeric_stitched_contigs += 1

        stats_dict[sample_name].append(str(chimeric_stitched_contigs))

        ################################################################################################################
        # Total bases recovered (not counting N characters):
        ################################################################################################################
        stats_dict[sample_name].append(str(sample_name_to_total_bases_dict[sample_name]))

    ####################################################################################################################
    # SeqLengths:
    ####################################################################################################################
    for sample_name in stats_dict:
        stats_dict_for_printing = '\t'.join(stats_dict[sample_name])
        lines_for_stats_report.append(f'{sample_name}\t{stats_dict_for_printing}')

    with open(f'{args.stats_filename}.tsv', 'w') as hybpiper_stats_handle:
        for item in lines_for_stats_report:
            if len([item for item in item.split('\t')]) == 2:  # i.e. no bam file and no stats
                continue
            hybpiper_stats_handle.write(f'{item}\n')

    logger.info(f'{"[INFO]:":10} A statistics table has been written to file: {args.stats_filename}.tsv')


if __name__ == "__main__":
    standalone()
