#!/usr/bin/env python

"""
usage: python distribute_targets.py baitfile\n

Given a file containing all of the "baits" for a target enrichment, create separate
FASTA files with all copies of that bait. Multiple copies of the same bait can be 
specified using a "-" delimiter. For example, the following will be sorted in to the same
file:

Anomodon-rbcl
Physcomitrella-rbcl

The results can come from either BLASTx or BWA.

Given multiple baits, the script will choose the most appropriate 'reference' sequence
using the highest cumulative BLAST scores or Mapping Quality across all hits.

Output directories can also be created, one for each target category (the default is to put them all in the current one)
The field delimiter may also be changed.
"""

import sys
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

    :param sequence:
    :return:
    """

    remainder = len(sequence) % 3
    return sequence if remainder == 0 else sequence + Seq('N' * (3 - remainder))


def mkdir_p(path):
    """

    :param path:
    :return:
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def tailored_target_blast(blastxfilename, exclude=None):
    """
    Determine, for each protein, the 'best' target protein, by tallying up the blastx hit scores.

    :param blastxfilename:
    :param exclude:
    :return:
    """
    blastxfile = open(blastxfilename)
    
    hitcounts = {}
    for result in blastxfile:
        result = result.split()
        hitname = result[1].split('-')
        bitscore = float(result[-1])
        try: 
            protname = hitname[1]
        except IndexError:
            raise IndexError('Gene name not found! FASTA headers should be formatted like this:\n >SpeciesName-GeneName\n')
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


def tailored_target_bwa(bamfilename, unpaired=False, exclude=None):
    """
    Determine, for each protein, the 'best' target protein, by tallying up the BWA map scores.

    :param bamfilename:
    :param unpaired:
    :param exclude:
    :return:
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
            raise IndexError('Gene name not found! FASTA headers should be formatted like this:\n >SpeciesName-GeneName\n')
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

     
def distribute_targets(baitfile, dirs, delim, besthits, translate=False, target=None):
    """

    :param baitfile:
    :param dirs:
    :param delim:
    :param besthits:
    :param translate:
    :param target:
    :return:
    """
    if target:
        if os.path.isfile(target):
            print(('[DISTRIBUTE]: Reading preferred target names from {} \n'.format(target)))
            genes_to_targets = {x.split()[0]: x.rstrip().split()[1] for x in open(target)}
            target_is_file = True
        else:
            target_is_file = False    
        
    targets = SeqIO.parse(baitfile, 'fasta')
    no_matches = []
    for prot in targets:
        # Get the 'basename' of the protein
        prot_cat = prot.id.split(delim)[-1]
        if translate:
            prot.seq = pad_seq(prot.seq).translate()  # seq padded to avoid BioPython warning
        
        if dirs:
            mkdir_p(prot_cat)
        
        if prot_cat in besthits:
            if target:
                if target_is_file:
                    besthit_taxon = genes_to_targets[prot_cat]
                else:
                    besthit_taxon = target
            else:       
                besthit_taxon = besthits[prot_cat]
            if prot.id.split("-")[0] == besthit_taxon:
                outfile = open(os.path.join(prot_cat, f'{prot_cat}_baits.fasta'), 'w')
                SeqIO.write(prot, outfile, 'fasta')
                outfile.close()
        else:
            no_matches.append(prot_cat)
    logger.info(f'[DISTRIBUTE]: {len(set(no_matches))} proteins had no good matches.')


def help():
    print(__doc__)
    return


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d', '--delimiter', help='Field separating FASTA ids for multiple sequences per target. '
                                                  'Default is '-' . For no delimeter, write None', default='-')
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
        besthits = tailored_target_blast(args.blastx, args.exclude)
        translate = False            
    if args.bam:
        translate = True
        besthits = tailored_target_bwa(args.bam, args.unpaired, args.exclude)

    distribute_targets(args.baitfile,
                       dirs=True,
                       delim=args.delimiter,
                       besthits=besthits,
                       translate=translate,
                       target=args.target)


if __name__ == '__main__':
    main()
