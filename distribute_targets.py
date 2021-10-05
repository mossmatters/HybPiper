#!/usr/bin/env python

"""
usage: python distribute_targets.py baitfile

Taks a file containing all of the "baits" for a target enrichment. The file can contain multiple copies of the same
bait as specified using a "-" delimiter. For example, the following:

Anomodon-rbcl
Physcomitrella-rbcl

Given multiple baits, the script will choose the most appropriate 'reference' sequence
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

# Create logger:
logger = logging.getLogger(f'__main__.{__name__}')


def pad_seq(sequence):
    """
    Pads a sequence Seq object to a multiple of 3 with 'N'.

    :param Bio.SeqRecord.SeqRecord sequence: sequence to pad
    :return: Bio.SeqRecord.SeqRecord sequence padded with Ns if required, bool True if sequence needed padding
    """

    remainder = len(sequence.seq) % 3
    if remainder == 0:
        return sequence, False
    else:
        # logger.info(f'{"[WARN!]:":10} The baitfile nucleotide sequence {sequence.id} is not a multiple of 3!')
        sequence.seq = sequence.seq + Seq('N' * (3 - remainder))
        return sequence, True


def mkdir_p(path):
    """
    Creates a directory corresponding the the given path, if it doesn't already exist.

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
    :param str exclude: no not use any target sequence specified by this string
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
        try: 
            protname = hitname[1]
        except IndexError:
            raise IndexError('Gene name not found! FASTA headers should be formatted like this:\n '
                             '>SpeciesName-GeneName\n')
        taxon = hitname[0]
        if exclude and exclude in taxon:
            continue
        else:
            if protname in hitcounts:
                if taxon in hitcounts[protname]:
                    hitcounts[protname][taxon] += bitscore
                else:
                    hitcounts[protname][taxon] = bitscore
            else:
                hitcounts[protname] = {taxon: 1}

    # For each protein, find the taxon with the highest total hit bitscore.
    besthits = {}
    besthit_counts = {}
    print(f'\nlength of hitcounts dict is: {len(hitcounts)}\n')
    for prot in hitcounts:
        top_taxon = sorted(iter(hitcounts[prot].keys()), key=lambda k: hitcounts[prot][k], reverse=True)[0]
        print(f'top taxon for prot {prot} is: \t{top_taxon}')
        besthits[prot] = top_taxon
        if top_taxon in besthit_counts:
            besthit_counts[top_taxon] += 1
        else:
            besthit_counts[top_taxon] = 1
    print(f'')
    tallyfile = open('bait_tallies.txt', 'w')
    print(f'length of besthit_counts dict is: {len(besthit_counts)}\n')
    for x in besthit_counts:
        tallyfile.write(f'{x}\t{besthit_counts[x]}\n')
    tallyfile.close()
    for key, value in besthits.items():
        print(key, value)
    print(f'\nlength of besthits dict is {len(besthits)}\n')
    return besthits        


def tailored_target_bwa(bamfilename, unpaired=False, exclude=None):
    """
    Determine, for each protein, the 'best' target protein, by tallying up the BWA map scores.

    :param str bamfilename: path to *.bam output of BWA mapping
    :param bool unpaired: if True, process a *_unpaired.bam file
    :param str exclude: no not use any target sequence specified by this string
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
        try: 
            protname = hitname[1]
        except IndexError:
            raise IndexError('Gene name not found! FASTA headers should be formatted like this:\n '
                             '>SpeciesName-GeneName\n')
        taxon = hitname[0]
        if exclude and exclude in taxon:
            continue
        else:
            if protname in hitcounts:
                if taxon in hitcounts[protname]:
                    hitcounts[protname][taxon] += mapscore
                else:
                    hitcounts[protname][taxon] = mapscore
            else:
                hitcounts[protname] = {taxon: 1}

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
    tallyfile = open('bait_tallies.txt', 'w')
    for x in besthit_counts:
        tallyfile.write(f'{x}\t{besthit_counts[x]}\n')
    tallyfile.close()
    return besthits    

     
def distribute_targets(baitfile, delim, besthits, translate=False, target=None):
    """
    Writes the single 'best' protein sequence from the target file (translated if neccessary) as a fasta file for each
    gene.

    :param str baitfile: path to baitfile
    :param str delim: symbol to use as gene delimeter; default is '-'
    :param dict besthits: dictionary of besthits[prot] = top_taxon
    :param bool translate: If True, translate nucleotide target SeqObject
    :param str target: always choose the target specified by this string
    :return:
    """

    if target:
        if os.path.isfile(target):
            logger.info(f'{"[NOTE]:":10} Reading preferred target names from {target}')
            genes_to_targets = {x.split()[0]: x.rstrip().split()[1] for x in open(target)}
            target_is_file = True
        else:
            target_is_file = False    
        
    targets = SeqIO.parse(baitfile, 'fasta')
    no_matches = []
    for sequence in targets:
        # Get the 'basename' of the protein
        gene_id = sequence.id.split(delim)[-1]
        if translate:
            seq, needed_padding = pad_seq(sequence)
            sequence.seq = seq.seq.translate()
        # print(prot)

        mkdir_p(gene_id)
        
        if gene_id in besthits:
            if target:
                if target_is_file:
                    besthit_taxon = genes_to_targets[gene_id]
                else:
                    besthit_taxon = target
            else:       
                besthit_taxon = besthits[gene_id]
            if sequence.id.split("-")[0] == besthit_taxon:
                with open(os.path.join(gene_id, f'{gene_id}_baits.fasta'), 'w') as ref_bait_seq_file:
                    SeqIO.write(sequence, ref_bait_seq_file, 'fasta')
        else:
            no_matches.append(gene_id)
    logger.info(f'{"[NOTE]:":10} {len(set(no_matches))} proteins had no good matches.')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d', '--delimiter', help='Field separating FASTA ids for multiple sequences per target. '
                                                  'Default is "-" . For no delimeter, write None', default='-')
    parser.add_argument('baitfile', help='FASTA file containing bait sequences')
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

    distribute_targets(args.baitfile,
                       delim=args.delimiter,
                       besthits=besthits,
                       translate=translate,
                       target=args.target)


if __name__ == '__main__':
    main()
