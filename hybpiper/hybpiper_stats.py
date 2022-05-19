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


def get_seq_lengths(targetfile, namelist, targetfile_sequence_type, sequence_type_to_calculate_stats_for,
                    seq_lengths_filename):
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
    :return str seq_lengths_report_filename: path to the sequence length report file written by this function
    """

    lines_for_report = []  # lines to write to file

    # Set variable 'filetype', used to reconstruct path to each sequence:
    if sequence_type_to_calculate_stats_for.upper() == 'GENE':
        filetype = 'FNA'
    elif sequence_type_to_calculate_stats_for.upper() == "SUPERCONTIG":
        filetype = 'supercontig'

    if not os.path.isfile(targetfile):
        print(f'Target file {targetfile} not found!')
        sys.exit()

    if not os.path.isfile(namelist):
        print(f'Name list file {namelist} not found!')
        sys.exit()

    namelist_parsed = [n.rstrip() for n in open(namelist).readlines()]

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

    # Get seq lengths for sample gene sequences (FNA, FAA or supercontigs):
    for name in namelist_parsed:
        parentDir, name = os.path.split(name)
        if not name:
            parentDir, name = os.path.split(parentDir)
        name_lengths = []  # lengths of sequences in nucleotides
        for gene in range(len(unique_names)):

            # Reconstruct path to the sequence:
            if filetype == 'supercontig':
                read_file = os.path.join(parentDir, name, unique_names[gene], name, 'sequences', 'intron',
                                         f'{unique_names[gene]}_supercontig.fasta')
            else:
                read_file = os.path.join(parentDir, name, unique_names[gene], name, "sequences", filetype,
                                         f'{unique_names[gene]}.{filetype}')

            if os.path.exists(read_file):
                seq_length = len(SeqIO.read(read_file, 'fasta').seq.ungap(gap='N'))
                if seq_length > 1.5 * avg_ref_lengths[gene] and filetype != 'supercontig':
                    sys.stderr.write(f'****WARNING! Sequence length for {name} is more than 50% longer than'
                                     f' {unique_names[gene]} reference!\n')
                name_lengths.append(str(seq_length))
            else:
                name_lengths.append("0")

        lengths_to_write = '\t'.join(name_lengths)
        lines_for_report.append(f'{name}\t{lengths_to_write}')

    # Write report file "seq_lengths.tsv"
    seq_lengths_report_filename = f'{seq_lengths_filename}.tsv'
    with open(seq_lengths_report_filename, 'w') as seq_lengths_handle:
        for item in lines_for_report:
            seq_lengths_handle.write(f'{item}\n')
    print(f'A sequence length table has been written to file: {seq_lengths_filename}.tsv')

    return seq_lengths_report_filename


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


def enrich_efficiency_blastx(blastxfilename, sample_name):
    """
    Parse BLASTX results to calculate enrichment efficiency

    :param str blastxfilename: path to the *.tsv BLASTx output filename for a given sample
    :param str sample_name: name of a given sample
    :return str, str, str: values for input reads, mapped reads, and percent mapped reads
    """

    reads_with_hits = [x.split()[0] for x in open(blastxfilename)]
    if os.path.isfile(blastxfilename.replace(".blastx", "_unpaired.blastx")):
        reads_with_hits += [x.split()[0] for x in open(blastxfilename.replace(".blastx", "_unpaired.blastx"))]
    mappedReads = len(set(reads_with_hits))

    # Recover total numbers of reads in input files:
    if os.path.exists(f'{sample_name}/total_input_reads_paired.txt'):
        with open(f'{sample_name}/total_input_reads_paired.txt', 'r') as paired_number:
            total_input_reads = int(paired_number.read().rstrip())
    elif os.path.exists(f'{sample_name}/total_input_reads_single.txt'):
        with open(f'{sample_name}/total_input_reads_single.txt', 'r') as single_number:
            total_input_reads = int(single_number.read().rstrip())
    else:
        raise ValueError(f'No file containing total input paired or single-end read count found!')

    if os.path.exists(f'{sample_name}/total_input_reads_unpaired.txt'):
        with open(f'{sample_name}/total_input_reads_unpaired.txt', 'r') as unpaired_number:
            total_input_reads = total_input_reads + int(unpaired_number.read().rstrip())

    try:
        pctMapped = 100 * mappedReads / total_input_reads
    except ZeroDivisionError:
        pctMapped = 0.0

    return str(total_input_reads), str(mappedReads), "{0:.1f}".format(pctMapped)


def enrich_efficiency_bwa(bamfilename):
    """
    Run and parse samtools flagstat output, return number of reads and number on target. Calculate percentage of
    reads mapped.

    :param bamfilename:
    :return str, str, str: values for input reads, mapped reads, and percent mapped reads:
    """

    samtools_cmd = f'samtools flagstat {bamfilename}'
    child = subprocess.Popen(samtools_cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    flagstat_results = [line for line in child.stdout.readlines()]
    for line in flagstat_results:
        if re.search('primary$', line):
            numReads = float(line.split()[0])
        if re.search(r'\bprimary mapped\b', line):
            mappedReads = float(line.split()[0])

    if os.path.isfile(bamfilename.replace(".bam", "_unpaired.bam")):
        unpaired_samtools_cmd = f'samtools flagstat {bamfilename.replace(".bam", "_unpaired.bam")}'
        unpaired_child = subprocess.Popen(unpaired_samtools_cmd, shell=True, stdout=subprocess.PIPE,
                                          universal_newlines=True)
        flagstat_results = [line for line in unpaired_child.stdout.readlines()]
        for line in flagstat_results:
            if re.search('primary$', line):
                numReads += float(line.split()[0])
            if re.search(r'\bprimary mapped\b', line):
                mappedReads += float(line.split()[0])
    try:
        pctMapped = 100 * mappedReads / numReads
    except ZeroDivisionError:
        pctMapped = 0.0

    return str(int(numReads)), str(int(mappedReads)), "{0:.1f}".format(pctMapped)


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

    # Set target file type and path:
    if args.targetfile_dna:
        targetfile = args.targetfile_dna
        targetfile_type = 'DNA'
    elif args.targetfile_aa:
        targetfile = args.targetfile_aa
        targetfile_type = 'protein'

    # Get sequence lengths for recovered genes, and write them to file:
    seq_lengths_file_path = get_seq_lengths(targetfile,
                                            args.namelist,
                                            targetfile_type,
                                            args.sequence_type,
                                            args.seq_lengths_filename)

    lines_for_stats_report = []

    categories = ["Name",
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
                  "GenesWithChimeraWarning"
                  ]

    categories_for_printing = '\t'.join(categories)
    lines_for_stats_report.append(categories_for_printing)

    seq_length_dict = seq_length_calc(seq_lengths_file_path)
    stats_dict = {}

    for line in open(args.namelist):
        name = line.rstrip()
        stats_dict[name] = []

        # Enrichment Efficiency
        bamfile = f'{name}/{name}.bam'
        blastxfile = f'{name}/{name}.blastx'
        if os.path.isfile(bamfile):
            stats_dict[name] += enrich_efficiency_bwa(bamfile)
        elif os.path.isfile(blastxfile):
            stats_dict[name] += enrich_efficiency_blastx(blastxfile, name)
        else:
            sys.stderr.write(f'No .bam or .blastx file found for {name}\n')

        # Recovery Efficiency
        stats_dict[name] += recovery_efficiency(name)
        stats_dict[name] += seq_length_dict[name]

        # Paralogs - long
        if os.path.isfile(f'{name}/{name}_genes_with_long_paralog_warnings.txt'):
            paralog_warns = file_len(f'{name}/{name}_genes_with_long_paralog_warnings.txt')
            stats_dict[name].append(str(paralog_warns))
        else:
            stats_dict[name].append("0")

        # Paralogs - by contig depth across query protein
        num_genes_paralog_warning_by_depth = 0
        if os.path.isfile(f'{name}/{name}_genes_with_paralog_warnings_by_contig_depth.csv'):
            with open(f'{name}/{name}_genes_with_paralog_warnings_by_contig_depth.csv') as paralogs_by_depth:
                lines = paralogs_by_depth.readlines()
                for gene_stats in lines:
                    stat = gene_stats.split(',')[3].strip()
                    if stat == 'True':
                        num_genes_paralog_warning_by_depth += 1
        stats_dict[name].append(str(num_genes_paralog_warning_by_depth))

        # Supercontig information:
        stitched_contig_produced = 0
        no_stitched_contig = 0
        # supercontig_no_trimming = 0
        # supercontig_with_trimming = 0
        stitched_contig_skipped = 0
        if os.path.isfile(f'{name}/{name}_genes_with_stitched_contig.csv'):
            with open(f'{name}/{name}_genes_with_stitched_contig.csv') as stitched_contig_stats:
                lines = stitched_contig_stats.readlines()
                for gene_stats in lines:
                    stat = gene_stats.split(',')[2]
                    # print(stat)
                    if re.search('single Exonerate hit', stat):
                        no_stitched_contig += 1
                    elif re.search('Stitched contig produced', stat):
                        stitched_contig_produced += 1
                    # elif re.search('NODE', stat):
                    #     supercontig_with_trimming += 1
                    # elif re.search('no contig trimming performed', stat):
                    #     supercontig_no_trimming += 1
                    elif re.search('Stitched contig step skipped', stat):
                        stitched_contig_skipped += 1
        # supercontigs_total = supercontig_no_trimming + supercontig_with_trimming
        stitched_contigs_produced_total = stitched_contig_produced
        stats_dict[name].append(str(no_stitched_contig))
        stats_dict[name].append(str(stitched_contigs_produced_total))
        # stats_dict[name].append(str(supercontig_with_trimming))
        stats_dict[name].append(str(stitched_contig_skipped))

        chimeric_stitched_contigs = 0
        if os.path.isfile(f'{name}/{name}_genes_derived_from_putative_chimeric_stitched_contig.csv'):
            with open(f'{name}/{name}_genes_derived_from_putative_chimeric_stitched_contig.csv') as \
                    chimeric_stitched_contig_stats:
                lines = chimeric_stitched_contig_stats.readlines()
                for gene_stats in lines:
                    stat = gene_stats.split(',')[2]
                    if re.search(' Chimera WARNING for stitched contig.', stat):
                        chimeric_stitched_contigs += 1
        stats_dict[name].append(str(chimeric_stitched_contigs))

    # SeqLengths
    for name in stats_dict:
        stats_dict_for_printing = '\t'.join(stats_dict[name])
        lines_for_stats_report.append(f'{name}\t{stats_dict_for_printing}')

    with open(f'{args.stats_filename}.tsv', 'w') as hybpiper_stats_handle:
        for item in lines_for_stats_report:
            hybpiper_stats_handle.write(f'{item}\n')
    print(f'A statistics table has been written to file: {args.stats_filename}.tsv')


if __name__ == "__main__":
    standalone()
