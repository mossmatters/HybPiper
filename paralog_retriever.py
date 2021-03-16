#!/usr/bin/env python

"""
This script will retrieve paralog nucleotide (CDS) sequences for a specified
gene in all samples located in namelist.txt. It writes all the (unaligned) sequences to stdout.
If a sample does not have paralogs for that gene, the sequence in the FNA directory is retrieved instead.
"""

import os, sys, argparse
from Bio import SeqIO
import logging
import fnmatch


########################################################################################################################
########################################################################################################################
# Configure logger:

# Create a custom logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Create handlers
# c_handler = logging.StreamHandler(sys.stdout)
existing_log_file_numbers = [int(file.split('_')[-1]) for file in os.listdir('.') if
                             fnmatch.fnmatch(file, '*.mylog*')]
if not existing_log_file_numbers:
    new_log_number = 1
else:
    new_log_number = sorted(existing_log_file_numbers)[-1] + 1
# f_handler = logging.FileHandler(f'logging_file.mylog_{new_log_number}', mode='w')

# Create formatters and add it to handlers
# c_format = logging.Formatter('%(message)s')
# f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# c_handler.setFormatter(c_format)
# f_handler.setFormatter(f_format)

# Add handlers to the logger
# logger.addHandler(c_handler)
# logger.addHandler(f_handler)


# f_handler = logging.FileHandler(f'logging_file.mylog_{new_log_number}', mode='w')
# f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# f_handler.setFormatter(f_format)
# logger.addHandler(f_handler)


########################################################################################################################
########################################################################################################################


def retrieve_seqs(path, name, gene):
    # print(f'CJJ path, name: {path}, {name}')
    genes_to_skip = []
    seqs_to_write = None
    try:
        with open(f'{name}/{name}_supercontigs_with_discordant_reads.csv') as discordant:
            lines = discordant.readlines()
            for line in lines:
                genes_to_skip.append(line.split(',')[1])
        logger.info(f'genes_to_skip from sample {name}: {genes_to_skip}')
    except FileNotFoundError:
        logger.info(f'No discordant reads summary file found for gene {gene} sample {name}!')

    ##################### Standard recovery: no genes skipped; goes to stderr ##########################################

    if os.path.isdir(os.path.join(path, name, gene, name, 'paralogs')):
        seqs_to_write = [x for x in
                         SeqIO.parse(os.path.join(path, name, gene, name, 'paralogs', '{}_paralogs.fasta'.format(gene)),
                                     'fasta')]
        num_seqs = str(len(seqs_to_write))
    elif os.path.isfile(os.path.join(path, name, gene, name, 'sequences', 'FNA', '{}.FNA'.format(gene))):
        seqs_to_write = SeqIO.read(os.path.join(path, name, gene, name, 'sequences', 'FNA', '{}.FNA'.format(gene)),
                                   'fasta')
        num_seqs = "1"

    if seqs_to_write:
        SeqIO.write(seqs_to_write, sys.stderr, 'fasta')
    else:
        num_seqs = "0"


    ##################### Modified recovery: genes flagged as hybrids are skipped; goes to stdout ######################

    if gene in genes_to_skip:
        logger.info(f'Skipping gene {gene} for sample {name}- potential hybrid sequence!')
        num_seqs = "0"
        return num_seqs
    else:
        if os.path.isdir(os.path.join(path, name, gene, name, 'paralogs')):
            seqs_to_write = [x for x in
                             SeqIO.parse(os.path.join(path, name, gene, name, 'paralogs', '{}_paralogs.fasta'.format(gene)),
                                         'fasta')]
            num_seqs = str(len(seqs_to_write))
        elif os.path.isfile(os.path.join(path, name, gene, name, 'sequences', 'FNA', '{}.FNA'.format(gene))):
            seqs_to_write = SeqIO.read(os.path.join(path, name, gene, name, 'sequences', 'FNA', '{}.FNA'.format(gene)),
                                       'fasta')
            num_seqs = "1"

        if seqs_to_write:
            SeqIO.write(seqs_to_write, sys.stdout, 'fasta')
        else:
            num_seqs = "0"
        return num_seqs


def main():
    parser = argparse.ArgumentParser(description=(__doc__), formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('namelist', help='Text file containing list of HybPiper output directories, one per line.')
    parser.add_argument('gene', help="Name of gene to extract paralogs")

    args = parser.parse_args()

    f_handler = logging.FileHandler(f'{args.gene}.mylog_{new_log_number}', mode='w')
    f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    f_handler.setFormatter(f_format)
    logger.addHandler(f_handler)


    namelist = [x.rstrip() for x in open(args.namelist)]
    num_seqs = []
    for name in namelist:
        path, name = os.path.split(name)
        if not name:
            path, name = os.path.split(path)
        num_seqs.append(retrieve_seqs(path, name, args.gene))
    # sys.stderr.write("{}\t{}\n".format(args.gene, "\t".join(num_seqs)))
    logger.info("{}\t{}\n".format(args.gene, "\t".join(num_seqs)))


if __name__ == "__main__": main()