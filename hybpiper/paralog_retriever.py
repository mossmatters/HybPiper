#!/usr/bin/env python

"""
This script will retrieve paralog nucleotide (CDS) sequences for a specified gene in all samples located in
namelist.txt. It writes the unaligned sequences to folders, which can be user specified; defaults are 'paralogs_all'
and 'paralogs_no_chimeras'.

If a sample does not have paralogs for that gene, the sequence in the FNA directory is retrieved instead.
"""

import os
import sys
import argparse
from Bio import SeqIO
import logging


# Create a custom logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def retrieve_seqs(path, name, gene, fasta_dir_all=None, fasta_dir_no_chimeras=None):
    chimeric_genes_to_skip = []
    seqs_to_write = None

    # Make output directories:
    if not os.path.isdir(fasta_dir_all):
        os.mkdir(fasta_dir_all)
    if not os.path.isdir(fasta_dir_no_chimeras):
        os.mkdir(fasta_dir_no_chimeras)

    try:
        with open(f'{name}/{name}_genes_derived_from_putative_chimera_supercontigs.csv') as chimeric:
            lines = chimeric.readlines()
            for line in lines:
                chimeric_genes_to_skip.append(line.split(',')[1])
        logger.info(f'genes_to_skip from sample {name}: {chimeric_genes_to_skip}')
    except FileNotFoundError:
        logger.info(f'No chimeric supercontig summary file found for gene {gene} sample {name}!')

    # Normal recovery of all sequences; writes to folder fasta_dir_all:
    if os.path.isdir(os.path.join(path, name, gene, name, 'paralogs')):
        seqs_to_write = [x for x in SeqIO.parse(os.path.join(path, name, gene, name, 'paralogs',
                                                             f'{gene}_paralogs.fasta'), 'fasta')]
        num_seqs = str(len(seqs_to_write))
    elif os.path.isfile(os.path.join(path, name, gene, name, 'sequences', 'FNA', f'{gene}.FNA')):
        seqs_to_write = SeqIO.read(os.path.join(path, name, gene, name, 'sequences', 'FNA', f'{gene}.FNA'), 'fasta')
        num_seqs = "1"

    if seqs_to_write:
        SeqIO.write(seqs_to_write, f'{fasta_dir_all}/{name}_{gene}_paralogs_all.fasta', 'fasta')
    else:
        num_seqs = "0"

    # Skip any putative chimeric supercontig sequences; writes to folder fasta_dir_no_chimeras:
    if gene in chimeric_genes_to_skip:
        logger.info(f'Skipping gene {gene} for sample {name} - putative chimeric supercontig sequence!')
        num_seqs = "0"
        return num_seqs
    else:
        if os.path.isdir(os.path.join(path, name, gene, name, 'paralogs')):
            seqs_to_write = [x for x in SeqIO.parse(os.path.join(path, name, gene, name, 'paralogs',
                                                                 f'{gene}_paralogs.fasta'), 'fasta')]
            num_seqs = str(len(seqs_to_write))
        elif os.path.isfile(os.path.join(path, name, gene, name, 'sequences', 'FNA', f'{gene}.FNA')):
            seqs_to_write = SeqIO.read(os.path.join(path, name, gene, name, 'sequences', 'FNA', f'{gene}.FNA'), 'fasta')
            num_seqs = "1"

        if seqs_to_write:
            SeqIO.write(seqs_to_write, f'{fasta_dir_no_chimeras}/{name}_{gene}_paralogs_nochimeras.fasta', 'fasta')
        else:
            num_seqs = "0"

    return num_seqs
            

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('namelist', help='Text file containing list of HybPiper output directories, one per line.')
    parser.add_argument('gene', help="Name of gene to extract paralogs")
    parser.add_argument('--fasta_dir_all', help='Specify directory for output FASTA files (ALL)',
                        default='paralogs_all')
    parser.add_argument('--fasta_dir_no_chimeras', help='Specify directory for output FASTA files (no putative '
                                                        'chimeric sequences)', default='paralogs_no_chimeras')
    
    args = parser.parse_args()
    
    namelist = [x.rstrip() for x in open(args.namelist)]
    num_seqs = []
    for name in namelist:
        path, name = os.path.split(name)
        if not name:
            path, name = os.path.split(path)
        num_seqs.append(retrieve_seqs(path, name, args.gene, args.fasta_dir_all, args.fasta_dir_no_chimeras))

    num_seqs_to_print = '\t'.join(num_seqs)
    sys.stderr.write(f'{args.gene}\t{num_seqs_to_print}\n')


if __name__ == "__main__":
    main()
