#!/usr/bin/env python
import sys,os,argparse
from Bio import SeqIO

helptext = '''Usage: 
    python retrieve_sequences.py targets.fasta sequence_dir aa/dna/intron/supercontig

This script will get the sequences generated from multiple runs of the HybSeqPipeline (reads_first.py).
Specify either a directory with all the HybPiper output directories or a file containing sequences of interest. 
It retreives all the gene names from the bait file used in the run of the pipeline.

You must specify whether you want the protein (aa) or nucleotide (dna) sequences.
You can also specify 'intron' to retreive the intron sequences, 
or 'supercontig' to get intron and exon sequences.

Will output unaligned fasta files, one per gene, to current directory.
'''

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("targetfile",help="FASTA File containing target sequences")
parser.add_argument("sample_names",help="Directory containing Hybpiper Output OR a file containing HybPiper output names, one per line")
parser.add_argument("sequence_type",help="Type of sequence to extract",choices=["dna","aa","intron","supercontig"])
parser.add_argument("--hybpiper_dir",help="Specify directory containing HybPiper output")
parser.add_argument("--fasta_dir",help="Specify directory for output FASTA files")

if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
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
    filename = 'supercontig'

#Use gene names parsed from a bait file.
baitfile = args.targetfile 
target_genes_dict = SeqIO.to_dict(SeqIO.parse(baitfile,'fasta'))
target_genes = list(set([x.id.split('-')[-1] for x in SeqIO.parse(baitfile,'fasta')]))

if os.path.isdir(args.sample_names):
    sampledir = args.sample_names
    sample_names = [x for x in os.listdir(sampledir) if os.path.isdir(os.path.join(sampledir,x)) and not x.startswith('.')]
else:
    sample_names = [x.rstrip() for x in open(args.sample_names)]
    if args.hybpiper_dir:
        sampledir = args.hybpiper_dir
    else:
        sampledir = '.'
if args.fasta_dir:
    fasta_dir = args.fasta_dir
else:
    fasta_dir = '.'


print("Retreiving {} genes from {} samples".format(len(target_genes),len(sample_names)))


for gene in target_genes:
    gene_seqs = []
    for rec in gene_seqs:
        rec.id = rec.id.split("-")[0]
        rec.description = ''
    for sample in sample_names:
        if seq_dir == 'intron':
            sample_path = os.path.join(sampledir,sample,gene,sample,'sequences',seq_dir,"{}_{}.fasta".format(gene,filename))
        else:
            sample_path = os.path.join(sampledir,sample,gene,sample,'sequences',seq_dir,gene+'.'+seq_dir)
        if os.path.isfile(sample_path):
            gene_seqs.append(SeqIO.read(sample_path,'fasta'))
    print("Found {} sequences for {}.".format(len(gene_seqs),gene))
    
    if seq_dir == 'intron':
        outfilename = "{}_{}.fasta".format(gene,filename)
    else:
        outfilename = gene + '.' + seq_dir
    
    SeqIO.write(gene_seqs,open(os.path.join(fasta_dir,outfilename),'w'),'fasta')