#!/usr/bin/env python

#########################
# HybPiper Stats Script #
#########################

"""
Writes a report file called "seq_lengths.tsv". The first line contains gene names. The second line contains the length
of the reference sequences (baits). If there are multiple baits per gene, the mean length is reported. All other rows
contain one sample per line.

Parses the "seq_lengths.tsv" and gathers additional statistics about the HybPiper run.

For an explanation of columns, see github.com/mossmatters/HybPiper/wiki
"""

import argparse
import os
import sys
import subprocess
import re
from Bio import SeqIO


def get_seq_lengths(baitfile, namelist, sequence_type):
    """

    :param str baitfile:
    :param str namelist:
    :param str sequence_type:
    :return str seq_lengths_report_filename: path to the sequence length report file written by this function
    """

    lines_for_report = []

    if sequence_type.upper() == 'DNA':
        filetype = 'FNA'
    elif sequence_type.upper() == 'AA':
        filetype = 'FAA'
    elif sequence_type.upper() == "SUPERCONTIG":
        filetype = "supercontig"
    else:
        print(__doc__)
        sys.exit()

    if not os.path.isfile(baitfile):
        print(f'Baitfile {baitfile} not found!')
        sys.exit()

    if not os.path.isfile(namelist):
        print(f'Name list file {namelist} not found!')
        sys.exit()

    namelist_parsed = [n.rstrip() for n in open(namelist).readlines()]

    # Get the names and lengths for each sequence in the bait/target file:
    gene_names = []
    reference_lengths = {}
    for prot in SeqIO.parse(baitfile, "fasta"):
        protname = prot.id.split("-")[-1]
        gene_names.append(protname)
        if protname in reference_lengths:
            reference_lengths[protname].append(len(prot.seq))
        else:
            reference_lengths[protname] = [len(prot.seq)]

    unique_names = list(set(gene_names))
    avg_ref_lengths = [(sum(reference_lengths[gene])/len(reference_lengths[gene])) for gene in unique_names]

    # Capture the unique gene names and average length for each gene to write to a report:
    unique_names_to_write = '\t'.join(unique_names)
    avg_ref_lengths_to_write = '\t'.join([str(x) for x in avg_ref_lengths])
    # sys.stdout.write(f'Species\t{unique_names_to_write}\nMeanLength\t{avg_ref_lengths_to_write}\n')
    lines_for_report.append(f'Species\t{unique_names_to_write}')
    lines_for_report.append(f'MeanLength\t{avg_ref_lengths_to_write}')

    # Get seq lengths for sample gene sequences (FNA, FAA or supercontigs):
    for name in namelist_parsed:
        parentDir, name = os.path.split(name)
        if not name:
            parentDir, name = os.path.split(parentDir)
        name_lengths = []
        for gene in range(len(unique_names)):
            if filetype == "supercontig":
                # read_file = os.path.join(parentDir, name, unique_names[gene], name, "sequences", "intron",
                #                          "{}_supercontig.fasta".format(unique_names[gene]))
                read_file = os.path.join(parentDir, name, unique_names[gene], name, 'sequences', 'intron',
                                         f'{unique_names[gene]}_supercontig_without_Ns.fasta')
            else:
                read_file = os.path.join(parentDir, name, unique_names[gene], name, "sequences", filetype,
                                         f'{unique_names[gene]}.{filetype}')

            if os.path.exists(read_file):

                if filetype == 'FNA' or filetype == 'supercontig':
                    # Strip any Ns inserted between gaps between Exonerate hits (with respect to the query protein):
                    seq_length = len(SeqIO.read(read_file, 'fasta').seq.ungap(gap='N'))
                elif filetype == 'FAA':
                    # Strip any Xs (translated Ns) inserted between gaps between Exonerate hits (with respect to the
                    # query protein):
                    seq_length = len(SeqIO.read(read_file, 'fasta').seq.ungap(gap='X'))

                if seq_length > 1.5 * avg_ref_lengths[gene] and filetype != "supercontig":
                    sys.stderr.write(f"****WARNING! Sequence length for {name} is more than 50% longer than"
                                     f" {unique_names[gene]} reference!\n")
                name_lengths.append(str(seq_length))
            else:
                name_lengths.append("0")

        lengths_to_write = '\t'.join(name_lengths)
        lines_for_report.append(f'{name}\t{lengths_to_write}')
        # sys.stdout.write(f'{name}\t{lengths_to_write}\n')

    # Write report file "seq_lengths.tsv"
    seq_lengths_report_filename = 'seq_lengths.tsv'
    with open(seq_lengths_report_filename, 'w') as seq_lengths_handle:
        for item in lines_for_report:
            seq_lengths_handle.write(f'{item}\n')

    return seq_lengths_report_filename


def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])


def enrich_efficiency_blastx(blastxfilename, name):
    """
    Parse BLASTX results to calculate enrichment efficiency
    """

    reads_with_hits = [x.split()[0] for x in open(blastxfilename)]
    if os.path.isfile(blastxfilename.replace(".blastx", "_unpaired.blastx")):
        reads_with_hits += [x.split()[0] for x in open(blastxfilename.replace(".blastx", "_unpaired.blastx"))]
    mappedReads = len(set(reads_with_hits))

    # Recover total numbers of reads in input files:
    if os.path.exists(f'{name}/total_input_reads_paired.txt'):
        with open(f'{name}/total_input_reads_paired.txt', 'r') as paired_number:
            total_input_reads = int(paired_number.read().rstrip())
    elif os.path.exists(f'{name}/total_input_reads_single.txt'):
        with open(f'{name}/total_input_reads_single.txt', 'r') as single_number:
            total_input_reads = int(single_number.read().rstrip())
    else:
        raise ValueError(f'No file containing total input paired or single-end read count found!')

    if os.path.exists(f'{name}/total_input_reads_unpaired.txt'):
        with open(f'{name}/total_input_reads_unpaired.txt', 'r') as unpaired_number:
            total_input_reads = total_input_reads + int(unpaired_number.read().rstrip())

    try:
        pctMapped = mappedReads / total_input_reads
    except ZeroDivisionError:
        pctMapped = 0.0
    return str(total_input_reads), str(mappedReads), "{0:.3f}".format(pctMapped)


def enrich_efficiency_bwa(bamfilename):
    """
    Run and parse samtools flagstat output, return number of reads and number on target
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
        pctMapped = mappedReads / numReads
    except ZeroDivisionError:
        pctMapped = 0.0
    return str(int(numReads)), str(int(mappedReads)), "{0:.3f}".format(pctMapped)


def recovery_efficiency(name):
    """
    Report the number of genes with mapping hits, contigs, and exon sequences
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
    From the output of get_seq_lengths.py, calculate the number of genes with seqs, and at least a pct of the
    reference length
    """

    seq_length_dict = {}
    with open(seq_lengths_fn) as seq_len:
        gene_names = seq_len.readline()
        # if blastx_adjustment:
        #     target_lengths = [float(value) * 3 for value in seq_len.readline().split()[1:]]
        # else:
        #     target_lengths = seq_len.readline().split()[1:]
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
    parser.add_argument("baitfile", help="FASTA file containing bait sequences for each gene. If there are multiple "
                                         "baits for a gene, the id must be of the form: >Taxon-geneName")
    parser.add_argument("namelist", help="text file with names of HybPiper output directories, one per line")
    parser.add_argument("sequence_type", help="Sequence type (dna, aa or supercontig) to recover lengths for",
                        choices=["dna", "DNA", "aa", "AA", "supercontig", "SUPERCONTIG"])
    # parser.add_argument("--blastx_adjustment", dest="blastx_adjustment", action='store_true',
    #                     help="Adjust stats for when blastx is used i.e. protein references, in cases where "
    #                          "get_seq_lengths.py has been run with parameter <dna> rather than <aa>", default=False)
    args = parser.parse_args()
    main(args)


def main(args):
    """
    Entry point for the assemble.py module.

    :param argparse.Namespace args:
    """

    # Get sequence lengths for recovered genes, and write them to file:
    seq_lengths_file_path = get_seq_lengths(args.baitfile, args.namelist, args.sequence_type)

    # Set blastx_adjustment boolean to adjust stats for when blastx is used i.e. protein references, in cases where
    # get_seq_lengths() has been run with parameter <dna> rather than <aa>:
    # if args.sequence_type

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
                  "Genesat150pct",
                  "ParalogWarningsLong",
                  "ParalogWarningsDepth",
                  "GenesWithoutSupercontigs",
                  "GenesWithSupercontigs",
                  "GenesWithSupercontigSkipped",
                  "GenesWithChimeraWarning"
                  ]

    categories_for_printing = '\t'.join(categories)
    lines_for_stats_report.append(categories_for_printing)
    # sys.stdout.write(f'{categories_for_printing}\n')

    # seq_length_dict = seq_length_calc(seq_lengths_file_path, args.blastx_adjustment)
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
        supercontig_produced = 0
        no_supercontig = 0
        # supercontig_no_trimming = 0
        # supercontig_with_trimming = 0
        supercontig_skipped = 0
        if os.path.isfile(f'{name}/{name}_genes_with_supercontigs.csv'):
            with open(f'{name}/{name}_genes_with_supercontigs.csv') as supercontig_stats:
                lines = supercontig_stats.readlines()
                for gene_stats in lines:
                    stat = gene_stats.split(',')[2]
                    # print(stat)
                    if re.search('single Exonerate hit', stat):
                        no_supercontig += 1
                    elif re.search('Supercontig produced', stat):
                        supercontig_produced += 1
                    # elif re.search('NODE', stat):
                    #     supercontig_with_trimming += 1
                    # elif re.search('no contig trimming performed', stat):
                    #     supercontig_no_trimming += 1
                    elif re.search('Supercontig step skipped', stat):
                        supercontig_skipped += 1
        # supercontigs_total = supercontig_no_trimming + supercontig_with_trimming
        supercontigs_total = supercontig_produced
        stats_dict[name].append(str(no_supercontig))
        stats_dict[name].append(str(supercontigs_total))
        # stats_dict[name].append(str(supercontig_with_trimming))
        stats_dict[name].append(str(supercontig_skipped))

        chimeric_supercontigs = 0
        if os.path.isfile(f'{name}/{name}_genes_derived_from_putative_chimera_supercontigs.csv'):
            with open(f'{name}/{name}_genes_derived_from_putative_chimera_supercontigs.csv') as \
                    chimeric_supercontig_stats:
                lines = chimeric_supercontig_stats.readlines()
                for gene_stats in lines:
                    stat = gene_stats.split(',')[2]
                    if re.search(' Chimera WARNING for supercontig.', stat):
                        chimeric_supercontigs += 1
        stats_dict[name].append(str(chimeric_supercontigs))

    # SeqLengths
    for name in stats_dict:
        stats_dict_for_printing = '\t'.join(stats_dict[name])
        lines_for_stats_report.append(f'{name}\t{stats_dict_for_printing}')
        # sys.stdout.write(f'{name}\t{stats_dict_for_printing}\n')

    with open('hybpiper_stats.tsv', 'w') as hybpiper_stats_handle:
        for item in lines_for_stats_report:
            hybpiper_stats_handle.write(f'{item}\n')


if __name__ == "__main__":
    standalone()
