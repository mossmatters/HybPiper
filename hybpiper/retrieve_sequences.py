#!/usr/bin/env python

"""
This script will get the sequences generated from multiple runs of the 'hybpiper assemble' command.
Specify either a directory with all the HybPiper output directories or a file containing sample names of interest.
It retrieves all the gene names from the bait file used in the run of the pipeline.

You must specify whether you want the protein (aa), nucleotide (dna) sequences.

You can also specify 'intron' to retrieve the intron sequences, or 'supercontig' to get intron and exon sequences.

Will output unaligned fasta files, one per gene.
"""

import os
import sys
import argparse
from Bio import SeqIO


def recover_sequences_from_all_samples(seq_dir, filename, target_genes, sample_names, hybpiper_dir=None,
                                       fasta_dir=None):
    """
    Recovers sequences (dna, amino acid, supercontig or intron) for all genes from all samples

    :param str seq_dir: directory to recover sequence from
    :param str filename: file name component used to reconstruct path to file
    :param list target_genes: list of unique gene names in the bait/target file
    :param str sample_names: directory of samples, or text file with list of sample names
    :param None or str hybpiper_dir: if provided, a path to the directory containing HybPiper output
    :param None or str fasta_dir: directory name for output files, default is current directory
    :return None:
    """

    if os.path.isdir(sample_names):
        sampledir = sample_names
        sample_names = [x for x in os.listdir(sampledir) if os.path.isdir(os.path.join(sampledir, x)) and not
        x.startswith('.')]
    else:
        sample_names = [x.rstrip() for x in open(sample_names)]
        if hybpiper_dir:
            sampledir = hybpiper_dir
        else:
            sampledir = '.'

    if fasta_dir:
        fasta_dir = fasta_dir
        if not os.path.isdir(fasta_dir):
            os.mkdir(fasta_dir)
    else:
        fasta_dir = '.'

    print(f'Retrieving {len(target_genes)} genes from {len(sample_names)} samples')
    for gene in target_genes:
        numSeqs = 0

        # Construct names for intron and supercontig output files:
        if seq_dir in ['intron', 'supercontig']:
            outfilename = f'{gene}_{filename}.fasta'
        else:
            outfilename = f'{gene}.{seq_dir}'

        with open(os.path.join(fasta_dir, outfilename), 'w') as outfile:
            for sample in sample_names:

                # Get path to the gene/intron/supercontig sequence:
                if seq_dir == 'intron':
                    sample_path = os.path.join(sampledir, sample, gene, sample, 'sequences', seq_dir,
                                               f'{gene}_{filename}.fasta')
                else:
                    sample_path = os.path.join(sampledir, sample, gene, sample, 'sequences', seq_dir,
                                               f'{gene}.{seq_dir}')
                try:
                    seq = next(SeqIO.parse(sample_path, 'fasta'))
                    SeqIO.write(seq, outfile, 'fasta')
                    numSeqs += 1
                # except FileNotFoundError or StopIteration:  # BioPython 1.80 returns StopIteration error?
                except FileNotFoundError:
                    pass
        print(f'Found {numSeqs} sequences for {gene}.')


def recover_sequences_from_one_sample(seq_dir,
                                      filename,
                                      target_genes,
                                      single_sample_name,
                                      fasta_dir=None):
    """
    Recovers sequences (dna, amino acid, supercontig or intron) for all genes from one sample

    :param str seq_dir: directory to recover sequence from
    :param str filename: file name component used to reconstruct path to file
    :param list target_genes: list of unique gene names in the bait/target file
    :param str single_sample_name: directory of a single sample
    :param None or str fasta_dir: directory name for output files, default is current directory
    :return None:
    """

    if not os.path.isdir(single_sample_name):
        sys.exit(f'Can not find a directory for sample {single_sample_name}, exiting...')

    # Create a user-supplied directory if provided, or write to the current directory if not:
    if fasta_dir:
        fasta_dir = fasta_dir
        if not os.path.isdir(fasta_dir):
            os.mkdir(fasta_dir)
    else:
        fasta_dir = '.'
    print(f'Retrieving {len(target_genes)} genes from sample {single_sample_name}...')
    sequences_to_write = []
    numSeqs = 0

    # Construct names for intron and supercontig output files:
    if seq_dir in ['intron', 'supercontig']:
        outfilename = f'{filename}.fasta'
    else:
        outfilename = f'{seq_dir}.fasta'
    for gene in target_genes:
        # Get path to the gene/intron/supercontig sequence:
        if seq_dir == 'intron':
            sample_path = os.path.join(single_sample_name, gene, single_sample_name, 'sequences', seq_dir,
                                       f'{gene}_{filename}.fasta')
        else:
            sample_path = os.path.join(single_sample_name, gene, single_sample_name, 'sequences', seq_dir,
                                       f'{gene}.{seq_dir}')
        try:
            seq = next(SeqIO.parse(sample_path, 'fasta'))
            seq.id = f'{seq.id}-{gene}'
            sequences_to_write.append(seq)
            numSeqs += 1
        # except FileNotFoundError or StopIteration:  # BioPython 1.80 returns StopIteration error?
        except FileNotFoundError:
            pass
    print(f'Found {numSeqs} sequences for sample {single_sample_name}.')

    with open(os.path.join(fasta_dir, outfilename), 'w') as outfile:
        SeqIO.write(sequences_to_write, outfile, 'fasta')


def standalone():
    """
    Used when this module is run as a stand-alone script. Parses command line arguments and runs function main().:
    """

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('targetfile', help='FASTA File containing target sequences')
    parser.add_argument('--sample_names',
                        help='Directory containing Hybpiper output OR a file containing HybPiper output names, '
                             'one per line', default=None)
    parser.add_argument('--single_sample_name',
                        help='A single sample name to recover sequences for', default=None)
    parser.add_argument('sequence_type', help='Type of sequence to extract', choices=['dna', 'aa', 'intron',
                                                                                      'supercontig'])
    parser.add_argument('--hybpiper_dir', help='Specify directory containing HybPiper output', default=None)
    parser.add_argument('--fasta_dir', help='Specify directory for output FASTA files', default=None)

    args = parser.parse_args()
    main(args)


def main(args):
    """
    Entry point for the assemble.py module.

    :param argparse.Namespace args:
    """

    # Check some command line parameters:
    if not args.sample_names and not args.single_sample_name:
        sys.exit(f'Please supply either the --sample_names or --single_sample_name flag and corresponding arguments!')
    if args.sample_names and args.single_sample_name:
        sys.exit(f'Please supply either the --sample_names or --single_sample_name flag and corresponding arguments!')

    # Set sequence directory name and file names:
    if args.sequence_type == 'dna':
        seq_dir = "FNA"
        filename = None
    elif args.sequence_type == 'aa':
        seq_dir = "FAA"
        filename = None
    elif args.sequence_type == 'intron':
        seq_dir = 'intron'
        filename = 'introns'
    elif args.sequence_type == 'supercontig':
        seq_dir = 'intron'
        filename = 'supercontig'

    # Use gene names parsed from a bait file.
    baitfile = args.targetfile
    target_genes = list(set([x.id.split('-')[-1] for x in SeqIO.parse(baitfile, 'fasta')]))

    # Recover sequences from all samples:
    if args.sample_names:
        recover_sequences_from_all_samples(seq_dir,
                                           filename,
                                           target_genes,
                                           args.sample_names,
                                           args.hybpiper_dir,
                                           args.fasta_dir)
    elif args.single_sample_name:
        recover_sequences_from_one_sample(seq_dir,
                                          filename,
                                          target_genes,
                                          args.single_sample_name,
                                          args.fasta_dir)


if __name__ == "__main__":
    standalone()
