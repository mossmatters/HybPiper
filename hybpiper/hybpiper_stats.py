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
                    namelist,
                    targetfile_sequence_type,
                    sequence_type_to_calculate_stats_for,
                    seq_lengths_filename,
                    compressed_sample_dict):
    """
    Recover the sequence length of each target file gene (calculated as mean length if a representative sequence from
    more than one taxon is provided for a given gene). Calculate the percentage length recovery for each gene,
    for each sample. If a protein target file was used, convert target gene lengths to the number of nucleotides
    (i.e. amino-acids x 3).

    :param str targetfile: path to the targetfile
    :param str namelist: path to the text file containing sample names
    :param str targetfile_sequence_type: sequence type in the target file ('DNA' or 'protein')
    :param str sequence_type_to_calculate_stats_for: gene (in nucleotides) or supercontig (in nucleotides)
    :param str seq_lengths_filename: optional filename for seq_lengths file. Default is seq_lengths.tsv
    :param dict compressed_sample_dict: dictionary containing *.tar.gz contents for compressed sample directories
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

    if not os.path.isfile(targetfile):
        logger.error(f'{"[ERROR]:":10} Target file {targetfile} not found!')
        sys.exit()

    if not os.path.isfile(namelist):
        logger.error(f'{"[ERROR]:":10} Name list file {namelist} not found!')
        sys.exit()

    namelist_parsed = [n.rstrip() for n in open(namelist).readlines() if n.rstrip()]

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
    sample_name_2_total_bases_dict = defaultdict(int)

    for name in namelist_parsed:  # iterate over samples

        # Check is the sample directory is a compressed tarball:
        compressed_sample = False
        if name in compressed_sample_dict:
            compressed_sample = True

        name_lengths = []  # lengths of sequences in nucleotides
        for gene in range(len(unique_names)):  # iterate over genes

            # Reconstruct path to the sequence:
            if filetype == 'supercontig':
                read_file = os.path.join(name, unique_names[gene], name, 'sequences', 'intron',
                                         f'{unique_names[gene]}_supercontig.fasta')
            else:
                read_file = os.path.join(name, unique_names[gene], name, "sequences", filetype,
                                         f'{unique_names[gene]}.{filetype}')

            # Get sequence details depending on compressed vs non-compressed:
            seq_length = None

            if compressed_sample:
                file_exists = True if read_file in compressed_sample_dict[name] else False
                size = compressed_sample_dict[name][read_file] if read_file in compressed_sample_dict[name] else 0
                if file_exists and size > 0:
                    seq_length = utils.get_compressed_seq_length(name,
                                                                 read_file)

            else:
                file_exists = True if os.path.isfile(read_file) else False
                size = os.path.getsize(read_file) if os.path.isfile(read_file) else 0
                if file_exists and size > 0:
                    seq_length = len(SeqIO.read(read_file, 'fasta').seq.replace('N', ''))

            if file_exists:
                if size == 0:
                    logger.warning(f'{"[WARNING]:":10} File {read_file} exists, but is empty! A length of zero will be '
                                   f'recorded for this sequence.\n')
                    name_lengths.append("0")
                else:
                    if seq_length > 1.5 * avg_ref_lengths[gene] and filetype != 'supercontig':
                        logger.warning(f'{"[WARNING]:":10} Sequence length for {name} is more than 50% longer than'
                                       f' {unique_names[gene]} reference!\n')
                    name_lengths.append(str(seq_length))
                    sample_name_2_total_bases_dict[name] += seq_length
            else:
                name_lengths.append("0")

        lengths_to_write = '\t'.join(name_lengths)
        lines_for_report.append(f'{name}\t{lengths_to_write}')

    # Write report file "seq_lengths.tsv"
    seq_lengths_report_filename = f'{seq_lengths_filename}.tsv'
    with open(seq_lengths_report_filename, 'w') as seq_lengths_handle:
        for item in lines_for_report:
            seq_lengths_handle.write(f'{item}\n')
    logger.info(f'{"[INFO]:":10} A sequence length table has been written to file: {seq_lengths_filename}.tsv')

    return seq_lengths_report_filename, sample_name_2_total_bases_dict


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
                             blastxfilename,
                             compressed_sample_bool=False,
                             blastxfile_unpaired_exists=False,
                             total_input_reads_paired_exists=False,
                             total_input_reads_single_exists=False,
                             total_input_reads_unpaired_exists=False):
    """
    Parse BLASTX results to calculate enrichment efficiency

    :param str sample_name: name of a given sample
    :param str blastxfilename: path to the *.tsv BLASTx output filename for a given sample
    :param bool compressed_sample_bool: True if sample directory is a compressed tarball
    :param bool blastxfile_unpaired_exists: True if an unpaired BLASTX exists for this sample
    :param bool total_input_reads_paired_exists:
    :param bool total_input_reads_single_exists:
    :param bool total_input_reads_unpaired_exists:
    :return str, str, str: values for input reads, mapped reads, and percent mapped reads
    """

    unpaired_blastx_file = blastxfilename.replace(".blastx", "_unpaired.blastx")
    total_input_reads_paired = f'{sample_name}/total_input_reads_paired.txt'
    total_input_reads_single = f'{sample_name}/total_input_reads_single.txt'
    total_input_reads_unpaired = f'{sample_name}/total_input_reads_unpaired.txt'

    # Process stats from the main BLAST file:
    if compressed_sample_bool:
        blastx_lines = utils.get_compressed_file_lines(sample_name,
                                                       blastxfilename)
    else:
        blastx_handle = open(blastxfilename, 'r')
        blastx_lines = blastx_handle.readlines()

    reads_with_hits = [x.split('\t')[0] for x in blastx_lines if x]

    if not compressed_sample_bool:
        blastx_handle.close()

    # If an unpaired BLASTX exists for this sample, process stats:
    if blastxfile_unpaired_exists:
        if compressed_sample_bool:
            blastx_unpaired_lines = utils.get_compressed_file_lines(sample_name,
                                                                    unpaired_blastx_file)
        else:
            blastx_unpaired_handle = open(unpaired_blastx_file, 'r')
            blastx_unpaired_lines = blastx_unpaired_handle.readlines()

        reads_with_hits += [x.split('\t')[0] for x in blastx_unpaired_lines if x]

        if not compressed_sample_bool:
            blastx_unpaired_handle.close()

    mapped_reads = len(set(reads_with_hits))

    # Recover total numbers of reads in input files:
    if total_input_reads_paired_exists:

        if compressed_sample_bool:
            total_input_reads_paired_lines = utils.get_compressed_file_lines(sample_name,
                                                                             total_input_reads_paired)
            total_input_reads = int(total_input_reads_paired_lines[0].rstrip())
        else:
            with open(total_input_reads_paired, 'r') as paired_number:
                total_input_reads = int(paired_number.read().rstrip())

    elif total_input_reads_single_exists:

        if compressed_sample_bool:
            total_input_reads_single_lines = utils.get_compressed_file_lines(sample_name,
                                                                             total_input_reads_single)
            total_input_reads = int(total_input_reads_single_lines[0].rstrip())
        else:
            with open(total_input_reads_single, 'r') as single_number:
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

    if blastxfile_unpaired_exists:

        if compressed_sample_bool:
            total_input_reads_unpaired_lines = utils.get_compressed_file_lines(sample_name,
                                                                               total_input_reads_unpaired)
            total_input_reads = total_input_reads + int(total_input_reads_unpaired_lines[0].rstrip())
        else:
            with open(total_input_reads_unpaired, 'r') as unpaired_number:
                total_input_reads = total_input_reads + int(unpaired_number.read().rstrip())

    try:
        pct_mapped = 100 * mapped_reads / total_input_reads
    except ZeroDivisionError:
        pct_mapped = 0.0

    return str(total_input_reads), str(mapped_reads), "{0:.1f}".format(pct_mapped)


def enrich_efficiency_bwa(name,
                          compressed_sample_bool=False,
                          bam_file_unpaired_exists=False):
    """
    Run and parse samtools flagstat output, return number of reads and number on target. Calculate percentage of
    reads mapped.

    :param str name: sample name
    :param bool compressed_sample_bool: True if sample directory is a compressed tarball
    :param bool bam_file_unpaired_exists: True if an unpaired bamfile exists for this sample
    :return str, str, str: values for input reads, mapped reads, and percent mapped reads:
    """

    bam_flagstats_tsv_file = f'{name}/{name}_bam_flagstat.tsv'
    unpaired_bam_flagstats_tsv_file = f'{name}/{name}_unpaired_bam_flagstat.tsv'
    num_reads = 0
    mapped_reads = 0

    # Process stats from the main bam file:
    if compressed_sample_bool:
        bam_flagstat_lines = utils.get_compressed_file_lines(name,
                                                             bam_flagstats_tsv_file)
    else:
        bam_flagstat_handle = open(bam_flagstats_tsv_file, 'r')
        bam_flagstat_lines = bam_flagstat_handle.readlines()

    for line in bam_flagstat_lines:
        if re.search('primary$', line):
            num_reads = float(line.split('\t')[0])
        if re.search(r'\bprimary mapped$', line):
            mapped_reads = float(line.split('\t')[0])

    if not compressed_sample_bool:
        bam_flagstat_handle.close()

    # If an unpaired bam exists for this sample, process stats:
    if bam_file_unpaired_exists:
        if compressed_sample_bool:
            bam_flagstat_unpaired_lines = utils.get_compressed_file_lines(name,
                                                                          unpaired_bam_flagstats_tsv_file)
        else:
            bam_flagstat_unpaired_handle = open(unpaired_bam_flagstats_tsv_file, 'r')
            bam_flagstat_unpaired_lines = bam_flagstat_unpaired_handle.readlines()

        for line in bam_flagstat_unpaired_lines:
            if re.search('primary$', line):
                num_reads += float(line.split('\t')[0])
            if re.search(r'\bprimary mapped\b', line):
                mapped_reads += float(line.split('\t')[0])

        if not compressed_sample_bool:
            bam_flagstat_unpaired_handle.close()

    try:
        pct_mapped = 100 * mapped_reads / num_reads
    except ZeroDivisionError:
        pct_mapped = 0.0

    return str(int(num_reads)), str(int(mapped_reads)), "{0:.1f}".format(pct_mapped)


def recovery_efficiency(name):
    """
    Reports the number of genes with mapping hits, contigs, and exon sequences

    :param str name: sample name
    :return list: a list containing the number of genes with contigs, exonerate hits, and assembled sequences
    """

    txt_files = ["spades_genelist.txt",
                 "exonerate_genelist.txt",
                 "genes_with_seqs.txt"
                 ]

    my_stats = []
    for txt in txt_files:
        if os.path.isfile(f'{name}/{txt}'):
            my_stats.append(file_len(f'{name}/{txt}'))
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
    Entry point for the assemble.py module.

    :param argparse.Namespace args:
    """

    logger.info(f'{"[INFO]:":10} HybPiper version {__version__} was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.info(f'{fill}\n')

    # Set target file type and path:
    targetfile = None
    targetfile_type = None

    if args.targetfile_dna:
        targetfile = args.targetfile_dna
        targetfile_type = 'DNA'
    elif args.targetfile_aa:
        targetfile = args.targetfile_aa
        targetfile_type = 'protein'

    assert targetfile
    assert targetfile_type

    # Parse namelist and parse any *.tar.gz samples:
    compressed_sample_dict = dict()

    for line in open(args.namelist):
        name = line.rstrip()
        if name:
            if re.search('/', name):
                sys.exit(f'{"[ERROR]:":10} A sample name must not contain '
                         f'forward slashes. The file {args.namelist} contains: {name}')

            compressed_sample = f'{name}.tar.gz'
            if os.path.isfile(compressed_sample):
                compressed_sample_dict[name] = utils.parse_compressed_sample(compressed_sample)

    # Get sequence lengths for recovered genes, and write them to file, along with total bases recovered:
    (seq_lengths_file_path,
     sample_name_2_total_bases_dict) = get_seq_lengths(targetfile,
                                                       args.namelist,
                                                       targetfile_type,
                                                       args.sequence_type,
                                                       args.seq_lengths_filename,
                                                       compressed_sample_dict)

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

    seq_length_dict = seq_length_calc(seq_lengths_file_path)
    stats_dict = {}

    for line in open(args.namelist):  # iterate over samples
        name = line.rstrip()
        if name:
            # print(f'name: {name}')
            stats_dict[name] = []

            # Check is the sample directory is a compressed tarball:
            compressed_sample_bool = False
            if name in compressed_sample_dict:
                compressed_sample_bool = True

            # Enrichment Efficiency
            bamfile = f'{name}/{name}.bam'
            bamfile_unpaired = f'{name}/{name}_unpaired.bam'
            blastxfile = f'{name}/{name}.blastx'
            blastxfile_unpaired = f'{name}/{name}_unpaired.blastx'
            total_input_reads_paired = f'{name}/total_input_reads_paired.txt'
            total_input_reads_single = f'{name}/total_input_reads_single.txt'
            total_input_reads_unpaired = f'{name}/total_input_reads_unpaired.txt'

            if compressed_sample_bool:
                bam_file_exists = True if bamfile in compressed_sample_dict[name] else False
                bam_file_unpaired_exists = True if bamfile_unpaired in compressed_sample_dict[name] else False
                blastxfile_exists = True if blastxfile in compressed_sample_dict[name] else False
                blastxfile_unpaired_exists = True if blastxfile_unpaired in compressed_sample_dict[name] else False
                total_input_reads_paired_exists = \
                    True if total_input_reads_paired in compressed_sample_dict[name] else False
                total_input_reads_single_exists = \
                    True if total_input_reads_single in compressed_sample_dict[name] else False
                total_input_reads_unpaired_exists = \
                    True if total_input_reads_unpaired in compressed_sample_dict[name] else False

                if not bam_file_exists and not blastxfile_exists:
                    fill = utils.fill_forward_slash(f'{"[WARNING]:":10} No *.bam or *.blastx file found for '
                                                    f'{name}. No statistics will be recovered. Please check the log '
                                                    f'file for this sample!',
                                                    width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                                    break_on_forward_slash=True)
                    logger.warning(f'{fill}')
                    continue

                if bam_file_exists:
                    stats_dict[name] += enrich_efficiency_bwa(name,
                                                              compressed_sample_bool=compressed_sample_bool,
                                                              bam_file_unpaired_exists=bam_file_unpaired_exists)
                else:
                    stats_dict[name] += enrich_efficiency_blastx(
                        name,
                        blastxfile,
                        compressed_sample_bool=compressed_sample_bool,
                        blastxfile_unpaired_exists=blastxfile_unpaired_exists,
                        total_input_reads_paired_exists=total_input_reads_paired_exists,
                        total_input_reads_single_exists=total_input_reads_single_exists,
                        total_input_reads_unpaired_exists=total_input_reads_unpaired_exists
                    )

            else:
                bam_file_exists = True if os.path.isfile(bamfile) else False
                bam_file_unpaired_exists = True if os.path.isfile(bamfile_unpaired) else False
                blastxfile_exists = True if os.path.isfile(blastxfile) else False
                blastxfile_unpaired_exists = True if os.path.isfile(blastxfile_unpaired) else False
                total_input_reads_paired_exists = True if os.path.isfile(total_input_reads_paired) else False
                total_input_reads_single_exists = True if os.path.isfile(total_input_reads_single) else False
                total_input_reads_unpaired_exists = True if os.path.isfile(total_input_reads_unpaired) else False

                if not bam_file_exists and not blastxfile_exists:
                    fill = utils.fill_forward_slash(f'{"[WARNING]:":10} No *.bam or *.blastx file found for '
                                                    f'{name}. No statistics will be recovered. Please check the log '
                                                    f'file for this sample!',
                                                    width=90, subsequent_indent=' ' * 11, break_long_words=False,
                                                    break_on_forward_slash=True)
                    logger.warning(f'{fill}')
                    continue

                if bam_file_exists:
                    stats_dict[name] += enrich_efficiency_bwa(name,
                                                              bam_file_unpaired_exists=bam_file_unpaired_exists)

                else:
                    stats_dict[name] += enrich_efficiency_blastx(
                        name,
                        blastxfile,
                        blastxfile_unpaired_exists=blastxfile_unpaired_exists,
                        total_input_reads_paired_exists=total_input_reads_paired_exists,
                        total_input_reads_single_exists=total_input_reads_single_exists,
                        total_input_reads_unpaired_exists=total_input_reads_unpaired_exists
                    )

            print(stats_dict[name])
    #
    #         # Recovery Efficiency
    #         stats_dict[name] += recovery_efficiency(name)
    #         stats_dict[name] += seq_length_dict[name]
    #
    #         # Paralogs - long
    #         if os.path.isfile(f'{name}/{name}_genes_with_long_paralog_warnings.txt'):
    #             paralog_warns = file_len(f'{name}/{name}_genes_with_long_paralog_warnings.txt')
    #             stats_dict[name].append(str(paralog_warns))
    #         else:
    #             stats_dict[name].append("0")
    #
    #         # Paralogs - by contig depth across query protein
    #         num_genes_paralog_warning_by_depth = 0
    #         if os.path.isfile(f'{name}/{name}_genes_with_paralog_warnings_by_contig_depth.csv'):
    #             with open(f'{name}/{name}_genes_with_paralog_warnings_by_contig_depth.csv') as paralogs_by_depth:
    #                 lines = paralogs_by_depth.readlines()
    #                 for gene_stats in lines:
    #                     stat = gene_stats.split(',')[3].strip()
    #                     if stat == 'True':
    #                         num_genes_paralog_warning_by_depth += 1
    #
    #         stats_dict[name].append(str(num_genes_paralog_warning_by_depth))
    #
    #         # Stitched contig information:
    #         stitched_contig_produced = 0
    #         no_stitched_contig = 0
    #         stitched_contig_skipped = 0
    #         if os.path.isfile(f'{name}/{name}_genes_with_stitched_contig.csv'):
    #             with open(f'{name}/{name}_genes_with_stitched_contig.csv') as stitched_contig_stats:
    #                 lines = stitched_contig_stats.readlines()
    #                 for gene_stats in lines:
    #                     stat = gene_stats.split(',')[2]
    #                     if re.search('single Exonerate hit', stat):
    #                         no_stitched_contig += 1
    #                     elif re.search('Stitched contig produced', stat):
    #                         stitched_contig_produced += 1
    #                     elif re.search('Stitched contig step skipped', stat):
    #                         stitched_contig_skipped += 1
    #
    #         stitched_contigs_produced_total = stitched_contig_produced
    #         stats_dict[name].append(str(no_stitched_contig))
    #         stats_dict[name].append(str(stitched_contigs_produced_total))
    #         stats_dict[name].append(str(stitched_contig_skipped))
    #
    #         chimeric_stitched_contigs = 0
    #         if os.path.isfile(f'{name}/{name}_genes_derived_from_putative_chimeric_stitched_contig.csv'):
    #             with open(f'{name}/{name}_genes_derived_from_putative_chimeric_stitched_contig.csv') as \
    #                     chimeric_stitched_contig_stats:
    #
    #                 lines = chimeric_stitched_contig_stats.readlines()
    #                 for gene_stats in lines:
    #                     stat = gene_stats.split(',')[2]
    #                     if re.search(' Chimera WARNING for stitched contig.', stat):
    #                         chimeric_stitched_contigs += 1
    #
    #         stats_dict[name].append(str(chimeric_stitched_contigs))
    #
    #         # Total bases recovered (not counting N characters):
    #         stats_dict[name].append(str(sample_name_2_total_bases_dict[name]))
    #
    # # SeqLengths
    # for name in stats_dict:
    #     stats_dict_for_printing = '\t'.join(stats_dict[name])
    #     lines_for_stats_report.append(f'{name}\t{stats_dict_for_printing}')
    #
    # with open(f'{args.stats_filename}.tsv', 'w') as hybpiper_stats_handle:
    #     for item in lines_for_stats_report:
    #         if len([item for item in item.split('\t')]) == 2:  # i.e. no bam file and no stats
    #             continue
    #         hybpiper_stats_handle.write(f'{item}\n')
    #
    # logger.info(f'{"[INFO]:":10} A statistics table has been written to file: {args.stats_filename}.tsv')


if __name__ == "__main__":
    standalone()
