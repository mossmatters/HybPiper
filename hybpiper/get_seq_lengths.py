#!/usr/bin/env python

"""
Usage: python baitfile.fasta namelist.txt dna/aa

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


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit()

    baitfile = sys.argv[1]
    namelistfile = sys.argv[2]
    sequenceType = sys.argv[3]

    if sequenceType.upper() == 'DNA':
        filetype = 'FNA'
    elif sequenceType.upper() == 'AA':
        filetype = 'FAA'
    elif sequenceType.upper() == "SUPERCONTIG":
        filetype = "supercontig"
    else:
        print(__doc__)
        sys.exit()

    if not os.path.isfile(baitfile):
        print(f'Baitfile {baitfile} not found!')
        sys.exit()

    if not os.path.isfile(namelistfile):
        print(f'Name list file {namelistfile} not found!')
        sys.exit()

    namelist = [n.rstrip() for n in open(namelistfile).readlines()]

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

    unique_names_to_write = '\t'.join(unique_names)
    avg_ref_lengths_to_write = '\t'.join([str(x) for x in avg_ref_lengths])
    sys.stdout.write(f'Species\t{unique_names_to_write}\nMeanLength\t{avg_ref_lengths_to_write}\n')

    for name in namelist:
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
                # Strip any Ns inserted between gaps between Exonerate hits (with respect to the query protein):
                seq_length = len(SeqIO.read(read_file, 'fasta').seq.ungap(gap='N'))

                if seq_length > 1.5 * avg_ref_lengths[gene] and filetype != "supercontig":
                    sys.stderr.write(f"****WARNING! Sequence length for {name} is more than 50% longer than"
                                     f" {unique_names[gene]} reference!\n")
                name_lengths.append(str(seq_length))
            else:
                name_lengths.append("0")

        lengths_to_write = '\t'.join(name_lengths)
        sys.stdout.write(f'{name}\t{lengths_to_write}\n')


if __name__ == '__main__':
    main()
