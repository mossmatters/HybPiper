#!/usr/bin/env python

#########################
# HybPiper Stats Script #
#########################

"""
Gather statistics about HybPiper run.

Supply the output of get_seq_lengths.py and a list of HybPiper directories

For an explanation of columns, see github.com/mossmatters/HybPiper/wiki
"""

import argparse
import os
import sys
import subprocess
import re


def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])


def enrich_efficiency_blastx(blastxfilename):
    """
    Parse BLASTX results to calculate enrichment efficiency
    """

    reads_with_hits = [x.split()[0] for x in open(blastxfilename)]
    if os.path.isfile(blastxfilename.replace(".blastx", "_unpaired.blastx")):
        reads_with_hits += [x.split()[0] for x in open(blastxfilename.replace(".blastx", "_unpaired.blastx"))]
    numreads = len(set(reads_with_hits))

    return "NA", str(numreads), "NA"  # TODO parse this info for the report


def enrich_efficiency_bwa(bamfilename):
    """
    Run and parse samtools flagstat output, return number of reads and number on target
    """

    samtools_cmd = f'samtools flagstat {bamfilename}'
    child = subprocess.Popen(samtools_cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    flagstat_results = [line for line in child.stdout.readlines()]
    numReads = float(flagstat_results[0].split()[0])
    mappedReads = float(flagstat_results[4].split()[0])

    if os.path.isfile(bamfilename.replace(".bam", "_unpaired.bam")):
        unpaired_samtools_cmd = f'samtools flagstat {bamfilename.replace(".bam", "_unpaired.bam")}'
        unpaired_child = subprocess.Popen(unpaired_samtools_cmd, shell=True, stdout=subprocess.PIPE,
                                          universal_newlines=True)
        flagstat_results = [line for line in unpaired_child.stdout.readlines()]
        numReads += float(flagstat_results[0].split()[0])
        mappedReads += float(flagstat_results[4].split()[0])
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


def seq_length_calc(seq_lengths_fn, blastx_adjustment):
    """
    From the output of get_seq_lengths.py, calculate the number of genes with seqs, and at least a pct of the
    reference length
    """

    seq_length_dict = {}
    with open(seq_lengths_fn) as seq_len:
        gene_names = seq_len.readline()
        if blastx_adjustment:
            target_lengths = [float(value) * 3 for value in seq_len.readline().split()[1:]]
        else:
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
    parser.add_argument("seq_lengths", help="Output of get_seq_lengths.py")
    parser.add_argument("namelist", help="Text file with names of HybPiper output directories, one per line")
    parser.add_argument("--blastx_adjustment", dest="blastx_adjustment", action='store_true',
                        help="Adjust stats for when blastx is used i.e. protein references, in cases where "
                             "get_seq_lengths.py has been run with parameter <dna> rather than <aa>", default=False)
    args = parser.parse_args()
    main(args)


def main(args):
    """
    Entry point for the hybpiper.py module.

    :param argparse.Namespace args:
    """

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
    sys.stdout.write(f'{categories_for_printing}\n')

    seq_length_dict = seq_length_calc(args.seq_lengths, args.blastx_adjustment)
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
            stats_dict[name] += enrich_efficiency_blastx(blastxfile)
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
        sys.stdout.write(f'{name}\t{stats_dict_for_printing}\n')


if __name__ == "__main__":
    standalone()
