#!/usr/bin/env python

"""
Usage:
    python retrieve_sequences.py targets.fasta sequence_dir aa/dna/intron/supercontig

This script will get the sequences generated from multiple runs of the HybSeqPipeline (reads_first.py).
Specify either a directory with all the HybPiper output directories or a file containing sequences of interest. 
It retreives all the gene names from the bait file used in the run of the pipeline.

You must specify whether you want the protein (aa) or nucleotide (dna) sequences.
You can also specify 'intron' to retreive the intron sequences, 
or 'supercontig' to get intron and exon sequences.

Will output unaligned fasta files, one per gene, to current directory.
"""

import os
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("targetfile", help="FASTA File containing target sequences")
parser.add_argument("sample_names", help="Directory containing Hybpiper Output OR a file containing HybPiper output "
                                         "names, one per line")
parser.add_argument("sequence_type", help="Type of sequence to extract", choices=["dna", "aa", "intron", "supercontig"])
parser.add_argument("--hybpiper_dir", help="Specify directory containing HybPiper output")
parser.add_argument("--fasta_dir", help="Specify directory for output FASTA files")

args = parser.parse_args()

if args.sequence_type == 'dna':
    seq_dir = "FNA"
elif args.sequence_type == 'aa':
    seq_dir = "FAA"
elif args.sequence_type == 'intron':
    seq_dir = 'intron'
    filename = 'introns'
elif args.sequence_type == 'supercontig':
    seq_dir = 'intron'
    # filename = 'supercontig'
    filename = 'intronerate_supercontig_without_Ns'

# Use gene names parsed from a bait file.
baitfile = args.targetfile
target_genes = list(set([x.id.split('-')[-1] for x in SeqIO.parse(baitfile, 'fasta')]))

if os.path.isdir(args.sample_names):
    sampledir = args.sample_names
    sample_names = [x for x in os.listdir(sampledir) if os.path.isdir(os.path.join(sampledir, x)) and not
                    x.startswith('.')]
else:
    sample_names = [x.rstrip() for x in open(args.sample_names)]
    if args.hybpiper_dir:
        sampledir = args.hybpiper_dir
    else:
        sampledir = '.'

if args.fasta_dir:
    fasta_dir = args.fasta_dir
    if not os.path.isdir(fasta_dir):
        os.mkdir(fasta_dir)
else:
    fasta_dir = '.'

print(f'Retrieving {len(target_genes)} genes from {len(sample_names)} samples')
for gene in target_genes:
    gene_seqs = []
    numSeqs = 0
    if seq_dir == 'intron':
        outfilename = f'{gene}_{filename}.fasta'
    else:
        outfilename = f'{gene}.{seq_dir}'

    with open(os.path.join(fasta_dir, outfilename), 'w') as outfile:
        for sample in sample_names:
            if seq_dir == 'intron':
                sample_path = os.path.join(sampledir, sample, gene, sample, 'sequences', seq_dir, f'{gene}'
                                                                                                  f'_{filename}.fasta')
                # print(f'intron sample path is {sample_path}')
            else:
                sample_path = os.path.join(sampledir, sample, gene, sample, 'sequences', seq_dir, f'{gene}.{seq_dir}')
            try:
                seq = next(SeqIO.parse(sample_path, 'fasta'))
                # print(f'seq is: {seq}')
                SeqIO.write(seq, outfile, 'fasta')
                numSeqs += 1
            # except FileNotFoundError or StopIteration:  # BioPython 1.80 returns StopIteration error?
            except FileNotFoundError:
                pass
    print(f'Found {numSeqs} sequences for {gene}.')
