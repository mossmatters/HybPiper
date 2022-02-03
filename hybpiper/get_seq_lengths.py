#!/usr/bin/env python

"""
Usage: python get-seq_lengths.py baitfile.fasta namelist.txt dna/aa

    Prepare a list of names from each of your samples that used the HybSeqPipeline, one name per line.
    Indicate whether the bait file is DNA or AA (affects length calculations).
    
    The script prints a table to stdout. The first line contains gene names.
    The second line contains the length of the reference sequences (baits). 
    If there are multiple baits per gene, the mean length is reported.
    All other rows contain one sample per line.
    
    This script requires BioPython to parse the bait FASTA file.
"""

import os
import sys
from Bio import SeqIO
import argparse


def standalone():
    """
    Used when this module is run as a stand-alone script. Parses command line arguments and runs function main().:
    """

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("baitfile", help="FASTA file containing bait sequences for each gene. If there are multiple "
                                         "baits for a gene, the id must be of the form: >Taxon-geneName")
    parser.add_argument("namelist", help="text file with names of HybPiper output directories, one per line")
    parser.add_argument("sequence_type", help="Sequence type (dna or aa) of the baitfile used")
    args = parser.parse_args()

    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit()

    main(args)


def main(args):
    """
    Entry point for the assemble.py module.

    :param argparse.Namespace args:
    """

    if args.sequence_type.upper() == 'DNA':
        filetype = 'FNA'
    elif args.sequence_type.upper() == 'AA':
        filetype = 'FAA'
    elif args.sequence_type.upper() == "SUPERCONTIG":
        filetype = "supercontig"
    else:
        print(__doc__)
        sys.exit()

    if not os.path.isfile(args.baitfile):
        print(f'Baitfile {args.baitfile} not found!')
        sys.exit()

    if not os.path.isfile(args.namelist):
        print(f'Name list file {args.namelist} not found!')
        sys.exit()

    namelist_parsed = [n.rstrip() for n in open(args.namelist).readlines()]

    # Get the names and lengths for each sequence in the bait/target file:
    gene_names = []
    reference_lengths = {}
    for prot in SeqIO.parse(args.baitfile, "fasta"):
        protname = prot.id.split("-")[-1]
        gene_names.append(protname)
        if protname in reference_lengths:
            reference_lengths[protname].append(len(prot.seq))
        else:
            reference_lengths[protname] = [len(prot.seq)]

    unique_names = list(set(gene_names))
    avg_ref_lengths = [(sum(reference_lengths[gene])/len(reference_lengths[gene])) for gene in unique_names]

    # Write the unique gene names and average length for each gene to stdout:
    unique_names_to_write = '\t'.join(unique_names)
    avg_ref_lengths_to_write = '\t'.join([str(x) for x in avg_ref_lengths])
    sys.stdout.write(f'Species\t{unique_names_to_write}\nMeanLength\t{avg_ref_lengths_to_write}\n')

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
                                         f'{unique_names[gene]}_intronerate_supercontig_without_Ns.fasta')
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
        sys.stdout.write(f'{name}\t{lengths_to_write}\n')


if __name__ == '__main__':
    standalone()
