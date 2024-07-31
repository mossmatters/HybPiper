#!/usr/bin/env python

"""
usage: python distribute_targets.py targetfile

Takes a file containing all the "targets" for a target enrichment. The file can contain multiple copies of the same
target as specified using a "-" delimiter. For example, the following:

Anomodon-rbcl
Physcomitrella-rbcl

Given multiple targets, the script will choose the most appropriate 'reference' sequence
using the highest cumulative BLAST scores or Mapping Quality (BWA) across all hits.

Output directories can also be created, one for each target category (the default is to put them all in the current one)
The field delimiter may also be changed.
"""

import os
import errno
import argparse
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import logging
from hybpiper import utils

# Create logger:
logger = logging.getLogger(f'hybpiper.assemble.{__name__}')


def mkdir_p(path):
    """
    Creates a directory corresponding to the given path, if it doesn't already exist.

    :param str path: path of directory to create
    :return:
    """

    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def tailored_target_blast(blastxfilename, unpaired=False, exclude=None):
    """
    Determine, for each protein, the 'best' target protein, by tallying up the blastx hit scores.

    :param str blastxfilename: path the BLASTx tabular output file
    :param bool unpaired: if True, process a *_unpaired.bam file
    :param str exclude: no not use any target sequence specified by this taxon name string
    :return: dict besthits: dictionary of besthits[prot] = top_taxon
    """

    with open(blastxfilename) as blastxfile_handle:
        blastx_results = blastxfile_handle.readlines()
    if unpaired:
        with open(blastxfilename.replace('.blastx', '_unpaired.blastx')) as blastxfile_unpaired_handle:
            blastx_results_unpaired = blastxfile_unpaired_handle.readlines()
        blastx_results += blastx_results_unpaired

    hitcounts = {}
    for result in blastx_results:
        result = result.split()
        hitname = result[1].split('-')
        bitscore = float(result[-1])
        protname = hitname[-1]
        taxon = '-'.join(hitname[:-1])
        if exclude and exclude in taxon:
            continue  # i.e. don't use sequences from this taxon
        else:
            if protname in hitcounts:
                if taxon in hitcounts[protname]:
                    hitcounts[protname][taxon] += bitscore
                else:
                    hitcounts[protname][taxon] = bitscore
            else:
                hitcounts[protname] = {taxon: bitscore}

    # For each protein, find the taxon with the highest total hit bitscore.
    besthits = {}
    besthit_counts = {}
    for prot in hitcounts:
        top_taxon = sorted(iter(hitcounts[prot].keys()), key=lambda k: hitcounts[prot][k], reverse=True)[0]
        besthits[prot] = top_taxon
        if top_taxon in besthit_counts:
            besthit_counts[top_taxon] += 1
        else:
            besthit_counts[top_taxon] = 1
    tallyfile = open('target_tallies.txt', 'w')
    for x in besthit_counts:
        tallyfile.write(f'{x}\t{besthit_counts[x]}\n')
    tallyfile.close()
    return besthits        


def tailored_target_bwa(bamfilename, unpaired=False, exclude=None):
    """
    Determine, for each protein, the 'best' target protein, by tallying up the BWA map scores.

    :param str bamfilename: path to *.bam output of BWA mapping
    :param bool unpaired: if True, process a *_unpaired.bam file
    :param str exclude: no not use any target sequence specified by this taxon name string
    :return: dict besthits: dictionary of besthits[prot] = top_taxon
    """

    samtools_cmd = 'samtools view -F 4 {}'.format(bamfilename)
    child = subprocess.Popen(samtools_cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    bwa_results = child.stdout.readlines()
    if unpaired:
        up_samtools_cmd = f'samtools view -F 4 {bamfilename.replace(".bam", "_unpaired.bam")}'
        up_child = subprocess.Popen(up_samtools_cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
        bwa_results += up_child.stdout.readlines()    
    hitcounts = {}
    for result in bwa_results:
        result = result.split()
        hitname = result[2].split('-')
        mapscore = float(result[4])
        protname = hitname[-1]
        taxon = '-'.join(hitname[:-1])

        if exclude and exclude in taxon:
            continue  # i.e. don't use sequences from this taxon
        else:
            if protname in hitcounts:
                if taxon in hitcounts[protname]:
                    hitcounts[protname][taxon] += mapscore
                else:
                    hitcounts[protname][taxon] = mapscore
            else:
                hitcounts[protname] = {taxon: mapscore}

    # For each protein, find the taxon with the highest total hit mapscore.
    besthits = {}
    besthit_counts = {}
    for prot in hitcounts:
        top_taxon = sorted(iter(hitcounts[prot].keys()), key=lambda k: hitcounts[prot][k], reverse=True)[0]
        besthits[prot] = top_taxon
        if top_taxon in besthit_counts:
            besthit_counts[top_taxon] += 1
        else:
            besthit_counts[top_taxon] = 1
    tallyfile = open('target_tallies.txt', 'w')
    for x in besthit_counts:
        tallyfile.write(f'{x}\t{besthit_counts[x]}\n')
    tallyfile.close()
    return besthits    

     
def distribute_targets(targetfile, delim, besthits, translate=False, target=None):
    """
    Writes the single 'best' protein sequence from the target file (translated if necessary) as a fasta file for each
    gene.

    :param str targetfile: path to targetfile
    :param str delim: symbol to use as gene delimeter; default is '-'
    :param dict besthits: dictionary of besthits[prot] = top_taxon
    :param bool translate: If True, translate nucleotide target SeqObject
    :param str target: always choose the target specified by this string (or target name in file)
    :return:
    """

    if target:
        if os.path.isfile(target):
            print(f'{"[NOTE]:":10} Reading preferred target names from {target}')
            logger.info(f'{"[NOTE]:":10} Reading preferred target names from {target}')
            genes_to_targets = {x.split()[0]: x.rstrip().split()[1] for x in open(target)}  # creates dictionary
            print(f'genes_to_targets is: {genes_to_targets}')
            target_is_file = True
        else:
            target_is_file = False    
        
    targets = SeqIO.parse(targetfile, 'fasta')
    no_matches = []
    for sequence in targets:
        # Get the 'basename' of the target sequence
        gene_id = sequence.id.split(delim)[-1]
        if translate:
            seq, needed_padding = utils.pad_seq(sequence)
            sequence.seq = seq.seq.translate()

        mkdir_p(gene_id)  # Make a directory for the gene
        
        if gene_id in besthits:
            if target:
                if target_is_file:
                    besthit_taxon = genes_to_targets[gene_id]  # recover specified target from dict using gene as key
                else:
                    besthit_taxon = target
            else:       
                besthit_taxon = besthits[gene_id]
            if '-'.join(sequence.id.split("-")[:-1]) == besthit_taxon:
                with open(os.path.join(gene_id, f'{gene_id}_target.fasta'), 'w') as ref_target_seq_file:
                    SeqIO.write(sequence, ref_target_seq_file, 'fasta')
        else:
            no_matches.append(gene_id)
    logger.info(f'{"[NOTE]:":10} {len(set(no_matches))} genes had no good matches.')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d', '--delimiter', help='Field separating FASTA ids for multiple sequences per target. '
                                                  'Default is "-" . For no delimeter, write None', default='-')
    parser.add_argument('targetfile', help='FASTA file containing target sequences')
    parser.add_argument('--blastx', help='tabular blastx results file, used to select the best target for each gene',
                        default=None)
    parser.add_argument('--bam', help='BAM file from BWA search, alternative to the BLASTx method', default=None)
    parser.add_argument('--target', help='Choose this version of the target always', default=None)
    parser.add_argument('--unpaired', help='Indicate whether to expect a file containing results from unpaired '
                                           'reads.', action='store_true', default=False)
    parser.add_argument('--exclude', help='Do not use any sequence with the specified string as the chosen target.',
                        default=None)
    args = parser.parse_args()

    if args.blastx:
        besthits = tailored_target_blast(args.blastx, args.unpaired, args.exclude)
        translate = False            
    if args.bam:
        translate = True
        besthits = tailored_target_bwa(args.bam, args.unpaired, args.exclude)

    distribute_targets(args.targetfile,
                       delim=args.delimiter,
                       besthits=besthits,
                       translate=translate,
                       target=args.target)


if __name__ == '__main__':
    main()
