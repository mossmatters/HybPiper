#!/usr/bin/env python

"""
This module takes a nucleotide SPAdes assembly for a gene (with contigs potentially containing exons and
introns/intergenic sequence) and a protein sequence for the 'best' target file sequence for the gene. It performs an
Exonerate search using the protein sequence as the query and the SPAdes contigs as targets (protein2genome model) and
extracts exon regions to create (where possible) a contiguous coding sequence for the gene/sample.

The module also performs the following tasks:

- Identification of introns and supercontigs (containing both exons and introns) using the function intronerate().
- Identification of paralogs.
- Identification of putative chimeric stitched contigs (derived from multiple paralogs).
"""

import shutil
import sys
import os
import subprocess
import argparse
import logging
from Bio import SeqIO, SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import tee
import re
from collections import defaultdict
from collections import Counter
from operator import itemgetter
import copy
import itertools
import shlex
from hybpiper import utils


def initial_exonerate(proteinfilename, assemblyfilename, prefix):
    """
    Conduct Exonerate search (first with option `--refine full` and then without if it doesn't work).

    Using the ryo option in Exonerate, the header contains useful information.

    https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-manual

    >%ti  = target id
    %qi   = query id
    %qab  = query alignment begin
    %qae  = query alignment end
    %pi   = percent id
    (%tS) = target strand
    %tab  = target alignment begin
    %tae  = target alignment end
    %tcs  = target coding sequence

    :param str proteinfilename: path to the chosen target-file protein query fasta file
    :param str assemblyfilename: path to the SPAdes assembly contigs file
    :param str prefix: path of gene/sample name
    :return NoneType or str, None or outputfilename:
        - None: returned if Exonerate searches fail to produce output file
        - outputfilename: the Exonerate text file output
    """

    logger = logging.getLogger(f'{os.path.split(prefix)[0]}')

    outputfilename = f'{prefix}/exonerate_results.fasta'
    exonerate_ryo = '">%ti,%qi,%qab,%qae,%pi,(%tS),%tab,%tae\\n%tcs\\n"'

    exonerate_command = f'exonerate -m protein2genome --showalignment yes --showvulgar no -V 0 --refine ' \
                        f'full --ryo' \
                        f' {exonerate_ryo} {proteinfilename} {assemblyfilename} > {outputfilename}'
    logger.debug(f'Exonerate command is: {exonerate_command}')

    try:  # with --refine
        subprocess.run(exonerate_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       universal_newlines=True, check=True)
    except subprocess.CalledProcessError as exc:
        logger.debug(f'Exonerate with "--refine" FAILED for {prefix}. Output is: {exc}')
        logger.debug(f'Exonerate with "--refine" stdout is: {exc.stdout}')
        logger.debug(f'Exonerate with "--refine" stderr is: {exc.stderr}')

    if utils.file_exists_and_not_empty(outputfilename):  # Exonerate with --refine can fail (empty output file) with no
        # error
        logger.debug('Exonerate ran with --refine')
        return outputfilename
    else:
        try:  # without --refine
            exonerate_command = f'exonerate -m protein2genome --showalignment yes --showvulgar no -V 0 ' \
                                f'--ryo' \
                                f' {exonerate_ryo} {proteinfilename} {assemblyfilename} > {outputfilename}'
            logger.debug(f'Exonerate command is: {exonerate_command}')
            subprocess.run(exonerate_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           universal_newlines=True, check=True)

        except subprocess.CalledProcessError as exc:
            logger.debug(f'Exonerate without "--refine" FAILED for {prefix}. Output is: {exc}')
            logger.debug(f'Exonerate without "--refine" stdout is: {exc.stdout}')
            logger.debug(f'Exonerate without "--refine" stderr is: {exc.stderr}')

    # Exonerate without --refine can fail (emtpy file) with no error:
    if utils.file_exists_and_not_empty(outputfilename):
        logger.debug('Exonerate ran without --refine')
        return outputfilename
    else:
        logger.debug('Exonerate failed without --refine; bad SPAdes contig(s)?')
        return None


def intronerate(exonerate_object, spades_contig_dict, logger=None, no_padding_supercontigs=False,
                keep_intermediate_files=False):
    """
    Attempts to identify introns within supercontigs, and writes fasta files containing 1) supercontig sequences
    containing exons AND introns, with 10 'N' characters inserted at any location SPAdes contigs have been
    concatenated; and 2) intron sequences only.

    :param __main__.Exonerate exonerate_object: object returned from the class Exonerate
    :param dict spades_contig_dict: a SeqIO dictionary of SPAdes contigs
    :param logging.Logger logger: a logger object
    :param bool no_padding_supercontigs: if True, don't pad contig joins in supercontigs with stretches if 10 Ns
    :param bool keep_intermediate_files: if True, keep intermediate files from supercontig
    :return NoneType: no explicit return value
    """

    logger.debug(f'no_padding_supercontigs: {no_padding_supercontigs}')

    # Make directory for Intronerate processing output:
    intronerate_processing_directory = f'{exonerate_object.prefix}/intronerate'
    if not os.path.exists(intronerate_processing_directory):
        os.mkdir(intronerate_processing_directory)

    # Make directory for writing Intronerate sequence output (intron and supercontig sequences):
    intronerate_sequence_directory = f'{exonerate_object.prefix}/sequences/intron'
    if not os.path.exists(intronerate_sequence_directory):
        os.mkdir(intronerate_sequence_directory)

    trimmed_hits_dict = exonerate_object.hits_subsumed_hits_removed_overlaps_trimmed_dict
    sample_name = os.path.split(exonerate_object.prefix)[-1]
    logger.debug(f'sample_name: {sample_name}')
    gene_name = os.path.split(exonerate_object.prefix)[-2]
    logger.debug(f'gene_name: {gene_name}')
    spades_contigs_for_intronerate_supercontig = []

    # Check whether there's more than one Exonerate hit for any given SPAdes contig:
    spades_contig_2_exonerate_hits_dict = defaultdict(list)
    all_exonerate_hit_contig_names_in_order = []
    for hit, hit_data_dict in trimmed_hits_dict.items():
        spades_name_only = hit.split(',')[0]
        all_exonerate_hit_contig_names_in_order.append(spades_name_only)
        spades_contig_2_exonerate_hits_dict[spades_name_only].append(hit_data_dict)

    count = Counter(all_exonerate_hit_contig_names_in_order)
    contigs_with_more_than_one_exonerate_hit = [key for key, value in count.items() if value >= 2]
    logger.debug(count)
    logger.debug(f'SPAdes contigs that have more than one Exonerate hit: {contigs_with_more_than_one_exonerate_hit}')

    # If there is more than one Exonerate hit for a given SPAdes contig, check that the hits are consecutive with
    # respect to the protein query, and that they all occur on the same strand; this is _expected_ to be the case,
    # or else something odd is going on:
    set_of_consecutive_spades_contigs = set()
    for i in range(len(trimmed_hits_dict) - 1):
        if all_exonerate_hit_contig_names_in_order[i] == all_exonerate_hit_contig_names_in_order[i + 1]:
            set_of_consecutive_spades_contigs.add(all_exonerate_hit_contig_names_in_order[i])

    logger.debug(f'Set of SPAdes contigs with Exonerate  hits that occur consecutively with respect to the protein '
                 f'query: {set_of_consecutive_spades_contigs}')

    if contigs_with_more_than_one_exonerate_hit:
        if set(contigs_with_more_than_one_exonerate_hit) != set_of_consecutive_spades_contigs:
            logger.info(f'\n{"[WARNING]:":10} There is more than one Exonerate hit for at least one SPAdes contig for '
                        f'gene {gene_name}, sample {sample_name}, but these do NOT appear consecutively with respect '
                        f'to the protein query. CHECK THIS! Intron and supercontig sequences will not be recovered.')
            return
        else:
            logger.debug(f'\nThere is more than one Exonerate hit for at least one SPAdes contig, and these appear '
                         f'consecutively with respect to the protein query. Proceeding...')

    # Check hit strands:
    for spades_contig in contigs_with_more_than_one_exonerate_hit:
        exonerate_dicts = spades_contig_2_exonerate_hits_dict[spades_contig]
        hit_strands = [exonerate_dict['hit_strand'] for exonerate_dict in exonerate_dicts]
        logger.debug(f'hit_strands for {spades_contig} are: {hit_strands}')
        if len(set(hit_strands)) != 1:
            logger.info(f'\n{"[WARNING]:":10} There is more than one Exonerate hit for {spades_contig} for gene'
                        f' {gene_name}, sample {sample_name}, but these do NOT appear on the same strand. CHECK THIS! '
                        f'Intron and supercontig sequences will not be recovered.')
            return
        else:
            logger.debug(f'There is more than one Exonerate hit for {spades_contig}, and these appear on the '
                         f'same strand. Proceeding...')

    # Process each hit:
    hit_spades_names = []  # list of contig names that have been added to the Intronerated supercontig list
    for hit, hit_data_dict in trimmed_hits_dict.items():
        spades_name_only = hit.split(',')[0]
        if spades_name_only in hit_spades_names:
            logger.debug(f'Contig {spades_name_only} has already been processed. Skipping this hit!')
            continue  # i.e. don't add this contig again
        else:
            hit_spades_names.append(spades_name_only)

        # Check if a hit corresponds to a SPAdes contig with more than one hit. If so, recover the hit with
        # the largest 3' trim coordinate. If 3' trimming was carried out, trim the corresponding SPAdes contig
        # according to this coordinate, and add the sequence to the Intronerated supercontig list:
        if spades_name_only in contigs_with_more_than_one_exonerate_hit:
            hit_with_largest_3prime_coordinate = max(spades_contig_2_exonerate_hits_dict[spades_name_only],
                                                     key=lambda hit_dict: hit_dict['query_range'][1])
            logger.debug(f'Contig {spades_name_only} has more than one Exonerate hit. Selecting the hit with the '
                         f'largest query_range coordinate: {hit_with_largest_3prime_coordinate}')
            hit_data_dict = hit_with_largest_3prime_coordinate

        # Get data for creating/writing reverse-complemented SeqRecord objects:
        hit_name = hit_data_dict['hit_sequence'].name
        hit_id = hit_data_dict['hit_sequence'].id
        hit_description = hit_data_dict['hit_sequence'].description

        hit_length = len(hit_data_dict['hit_sequence'])
        raw_spades_contig = spades_contig_dict[spades_name_only]
        raw_spades_contig_length = len(raw_spades_contig)
        raw_spades_contig_id = raw_spades_contig.id

        trimmed_hit_ranges_all = hit_data_dict['hit_range_all']
        inter_ranges_all = hit_data_dict['hit_inter_ranges']
        # print(hit_data_dict['hit_sequence'].description)
        three_prime_bases_trimmed = hit_data_dict['hit_sequence'].description.split(':')[-1].strip()

        # If no trimming has been performed for a SPAdes contig, or the overlap between adjacent Exonerate contig
        # hits is <= 3 amino acids, add the whole contig:
        if three_prime_bases_trimmed == 'N/A' or int(three_prime_bases_trimmed) <= 9:
            logger.debug(f'No trimming performed for {hit}, or trimming is <= 9 bases. Adding whole SPAdes contig to '
                         f'list! three_prime_bases_trimmed is: {three_prime_bases_trimmed}')
            if hit_data_dict['hit_strand'] == -1:
                revcomp_seqrecord = raw_spades_contig.reverse_complement()
                revcomp_seqrecord.name = hit_name
                revcomp_seqrecord.id = hit_id
                revcomp_seqrecord.description = hit_description
                spades_contigs_for_intronerate_supercontig.append(revcomp_seqrecord)
                continue
            else:
                spades_contigs_for_intronerate_supercontig.append(raw_spades_contig)
                continue

        # If trimming WAS performed, check if it's a multi-exon hit with intron(s) or frameshifts(?) present:
        if not inter_ranges_all:
            multi_exon_hit = False
        else:
            multi_exon_hit = True

        if not multi_exon_hit:
            logger.debug(f'Hit {hit} is a SINGLE-exon sequence! Trimming: {three_prime_bases_trimmed}, '
                         f'ranges: {trimmed_hit_ranges_all}')

            if hit_data_dict['hit_strand'] == 1:
                slice_coordinate = trimmed_hit_ranges_all[0][1] - int(three_prime_bases_trimmed)
                spades_contigs_for_intronerate_supercontig.append(raw_spades_contig[:slice_coordinate])

            elif hit_data_dict['hit_strand'] == -1:
                trimmed_hit_ranges_all = Exonerate.convert_coords_revcomp(trimmed_hit_ranges_all,
                                                                          raw_spades_contig_length)
                raw_spades_contig = raw_spades_contig.reverse_complement()
                raw_spades_contig.id = raw_spades_contig_id
                raw_spades_contig.name = raw_spades_contig_id
                raw_spades_contig.description = raw_spades_contig_id
                slice_coordinate = trimmed_hit_ranges_all[0][1] - int(three_prime_bases_trimmed)
                spades_contigs_for_intronerate_supercontig.append(raw_spades_contig[:slice_coordinate])

        elif multi_exon_hit:
            logger.debug(f'multi_exon_hit: {multi_exon_hit}')
            slice_found = False

            if hit_data_dict['hit_strand'] == -1:
                trimmed_hit_ranges_all = Exonerate.convert_coords_revcomp(trimmed_hit_ranges_all,
                                                                          raw_spades_contig_length)
                raw_spades_contig = raw_spades_contig.reverse_complement()
                raw_spades_contig.id = raw_spades_contig_id
                raw_spades_contig.name = raw_spades_contig_id
                raw_spades_contig.description = raw_spades_contig_id

            # Note that where Exonerate hits have been trimmed at their 3' end due to overlaps, the corresponding
            #  SPAdes contig can't simply be trimmed by the same number of bases, as it can comprise both exons and
            #  introns. Therefore, it's necessary to keep track of the cumulative length of the hit sequence within
            #  exons in each contig, and trim the contig based on the coordinate of the 3' end of the final exon in the
            #  Exonerate hit.
            cumulative_hit_span = 0
            # Check location of 3' trimmed hit end in hit ranges:
            for hit_range in trimmed_hit_ranges_all:
                logger.debug(f'hit_range: {hit_range}')
                # Get value for hit_length_adjusted; needs adjusting for hit location within a contig containing
                # introns as well as exons:
                hit_length_adjusted = hit_range[0] + (hit_length - cumulative_hit_span)
                logger.debug(f'Hit_length_adjusted: {hit_length_adjusted}')
                hit_span = hit_range[1] - hit_range[0]
                if hit_length_adjusted not in range(hit_range[0], hit_range[1]):  # i.e. 3' end of hit not in range
                    logger.debug(f'hit_length_adjusted {hit_length_adjusted} is NOT in range {hit_range[0]} -'
                                 f' {hit_range[1]}; the 3prime end of the stitched contig does not occur in this exon!')
                    cumulative_hit_span += hit_span  # keep track of stitched contig length covered by previous ranges
                else:
                    logger.debug(f'hit_length_adjusted {hit_length_adjusted} is FOUND in range {hit_range[0]} -'
                                 f' {hit_range[1]}. The 3prime end of the stitched contig occurs in this exon!')
                    slice_coordinate = hit_length_adjusted
                    spades_contigs_for_intronerate_supercontig.append(raw_spades_contig[:slice_coordinate])
                    slice_found = True
                    break
            if not slice_found:
                raise ValueError(f'Slice coordinate not found for hit {hit}')

    # Write concatenated Intronerated supercontig_seqrecord with Ns between hits:
    intron_supercontig_id = f'{sample_name}'

    intronerate_supercontig_seq_with_n = Seq('NNNNNNNNNN'.join([str(seq.seq) for seq in
                                                                spades_contigs_for_intronerate_supercontig]))
    intronerated_supercontig_seqrecord = SeqRecord(seq=intronerate_supercontig_seq_with_n,
                                                   id=intron_supercontig_id,
                                                   description='Joins between unique SPAdes contigs are separated by '
                                                               '10 "N" characters')
    with open(f'{intronerate_processing_directory}/{gene_name}_supercontig.fasta',
              'w') as intronerate_supercontig_handle:
        SeqIO.write(intronerated_supercontig_seqrecord, intronerate_supercontig_handle, 'fasta')

    # Write concatenated Intronerated supercontig_seqrecord without Ns between hits:
    intronerate_supercontig_seq = Seq(''.join([str(seq.seq) for seq in
                                               spades_contigs_for_intronerate_supercontig]))
    intronerated_supercontig_seqrecord = SeqRecord(seq=intronerate_supercontig_seq,
                                                   id=f'{intron_supercontig_id}_without_Ns',
                                                   description='Joins between unique SPAdes contigs are NOT separated '
                                                               'by 10 "N" characters')
    with open(f'{intronerate_processing_directory}/{gene_name}_supercontig_without_Ns.fasta',
              'w') as intronerate_supercontig_handle:
        SeqIO.write(intronerated_supercontig_seqrecord, intronerate_supercontig_handle, 'fasta')

    # Write individual constituent SPAdes contigs:
    with open(f'{intronerate_processing_directory}/{gene_name}_intronerate_supercontig_individual_contig_hits.fasta',
              'w') as intronerate_supercontig_hits_handle:
        SeqIO.write(spades_contigs_for_intronerate_supercontig, intronerate_supercontig_hits_handle, 'fasta')

    # Run Exonerate to get the gff file:
    intronerate_query = f'{exonerate_object.prefix}/sequences/FAA/{gene_name}.FAA'

    # Strip any "X" characters (inserted if sequence missing relative to target file reference):
    intronerate_query_stripped = f'{intronerate_processing_directory}/intronerate_query_stripped.fasta'
    with open(intronerate_query, 'r') as query_handle:
        query = SeqIO.read(query_handle, 'fasta')
        query.seq = query.seq.ungap('X')
        with open(intronerate_query_stripped, 'w') as stripped_handle:
            SeqIO.write(query, stripped_handle, 'fasta')

    if no_padding_supercontigs:
        ref = f'{gene_name}_supercontig_without_Ns.fasta'
        logger.debug(f'Intronerate run of Exonerate for gff will be run using {ref}')
        exonerate_supercontig_reference = f'{intronerate_processing_directory}/{ref}'
    else:
        ref = f'{gene_name}_supercontig.fasta'
        logger.debug(f'Intronerate run of Exonerate for gff will be run using {ref}')
        exonerate_supercontig_reference = f'{intronerate_processing_directory}/{ref}'

    exonerate_command = f'exonerate -m protein2genome -q {intronerate_query_stripped} -t ' \
                        f'{exonerate_supercontig_reference} ' \
                        f'--verbose 0 --showalignment yes --showvulgar no --refine full --showtargetgff yes >' \
                        f' {intronerate_processing_directory}/{gene_name}_intronerate_fasta_and_gff.txt'

    logger.debug(f'Intronerate run of Exonerate for gff with --refine: {exonerate_command}')

    try:  # with refine
        result = subprocess.run(exonerate_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True, check=True)
        logger.debug(f'Exonerate with "--refine" check_returncode() is: {result.check_returncode()}')
        logger.debug(f'Exonerate with "--refine" stdout is: {result.stdout}')
        logger.debug(f'Exonerate with "--refine" stderr is: {result.stderr}')
    except subprocess.CalledProcessError as exc:
        logger.debug(f'Exonerate with "--refine" FAILED for {gene_name}. Output is: {exc}')
        logger.debug(f'Exonerate with "--refine" stdout is: {exc.stdout}')
        logger.debug(f'Exonerate with "--refine" stderr is: {exc.stderr}')

    # Exonerate with --refine can fail (empty output file) with no error
    if utils.file_exists_and_not_empty(f'{intronerate_processing_directory}/{gene_name}_intronerate_fasta_and_gff.txt'):
        logger.debug('Exonerate run for Intronerate ran with --refine')
    else:
        exonerate_command = f'exonerate -m protein2genome -q {intronerate_query_stripped} -t ' \
                            f'{exonerate_supercontig_reference} --verbose 0 --showalignment yes --showvulgar no ' \
                            f'--showtargetgff yes > {intronerate_processing_directory}' \
                            f'/{gene_name}_intronerate_fasta_and_gff.txt'
        logger.debug(f'Intronerate run of Exonerate for gff without --refine: {exonerate_command}')
        try:  # without refine
            result = subprocess.run(exonerate_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True, check=True)
            logger.debug(f'Exonerate without "--refine" check_returncode() is: {result.check_returncode()}')
            logger.debug(f'Exonerate without "--refine" stdout is: {result.stdout}')
            logger.debug(f'Exonerate without "--refine" stderr is: {result.stderr}')
        except subprocess.CalledProcessError as exc:
            logger.debug(f'Exonerate without "--refine" FAILED for {gene_name}. Output is: {exc}')
            logger.debug(f'Exonerate without "--refine" stdout is: {exc.stdout}')
            logger.debug(f'Exonerate without "--refine" stderr is: {exc.stderr}')

    # Exonerate without --refine can fail (emtpy file) with no error
    if utils.file_exists_and_not_empty(f'{intronerate_processing_directory}/{gene_name}_intronerate_fasta_and_gff.txt'):
        logger.debug('Exonerate ran without --refine')
    else:
        logger.debug(f'Exonerate within intronerate() FAILED - skipping Intronerate for gene {gene_name}.')
        return

    # CJJ only use introns that occur in between detected exons, as there's no guarantee that the query protein is
    #  full length (we don't want flanking exonic regions in the raw SPAdes contigs being annotated as introns by
    #  mistake).

    # Parse out the gff lines only from the `intronerate_fasta_and_gff.gff` file:
    with open(f'{intronerate_processing_directory}/{gene_name}_intronerate_fasta_and_gff.txt',
              'r') as intronerate_fasta_gff_handle:
        gff_dump = intronerate_fasta_gff_handle.read()
        gff_split = re.split('# --- END OF GFF DUMP ---|# --- START OF GFF DUMP ---\n#\n#\n', gff_dump)
        gff_data_only = gff_split[1]  # TODO check if first hit result is always the best!
    with open(f'{intronerate_processing_directory}/intronerate.gff', 'w') as intronerate_gff_handle:
        intronerate_gff_handle.write(gff_data_only)

    # Parse the C4 alignment in file `intronerate_fasta_and_gff.txt` and write 1) Intronerated supercontigs; 2) introns
    # only:
    intronerate_supercontig = SeqIO.read(exonerate_supercontig_reference, 'fasta')

    exonerate_searchio_alignment = list(SearchIO.parse(
        f'{intronerate_processing_directory}/{gene_name}_intronerate_fasta_and_gff.txt', 'exonerate-text'))

    single_exonerate_qresult = exonerate_searchio_alignment[0]

    if len(single_exonerate_qresult.hsps) != 1:
        logger.debug(f'searchio_object hsps list is greater than 1 for {gene_name}!')
        # raise ValueError(f'searchio_object hsps list is greater than 1!')
    single_hsp = single_exonerate_qresult.hsps[0]  # TODO fix issue e.g. DEA_14928_S84 gene 6978 has two real HSPs
    intron_sequences = []
    for inter_range in single_hsp.hit_inter_ranges:
        intron_seqrecord = intronerate_supercontig[inter_range[0]:inter_range[1]]
        intron_seqrecord.description = 'intron'
        intron_sequences.append(intron_seqrecord)

    # Only write an "introns.fasta" file in certain cases:
    if exonerate_object.stitched_contig_seqrecord.description == 'single_hit' and \
            len(exonerate_object.hits_subsumed_hits_removed_overlaps_trimmed_dict['hit_inter_ranges']) == 0:
        logger.debug(f'Sequence for gene {gene_name} is derived from a single Exonerate hit with no introns - '
                     f'only a supercontig sequence will be recovered (i.e. no "introns.fasta" file) will be recovered')
    elif len(intron_sequences) == 0:  # e.g. When Exonerate splits the intronerate_supercontig target
        logger.debug(f'No introns for gene {gene_name}! This is likely caused by Exonerate splitting the '
                     f'"intronerate_supercontig" target in to multiple HSPs (i.e. it could not find an '
                     f'intron). Skipping intron recovery for this gene.')
    else:
        with open(f'{intronerate_processing_directory}/{gene_name}_introns.fasta', 'w') as introns_fasta_handle:
            SeqIO.write(intron_sequences, introns_fasta_handle, 'fasta')

        # Move the intron sequence to the intronerate_sequence_directory:
        if os.path.exists(f'{intronerate_sequence_directory}/{gene_name}_introns.fasta'):
            os.remove(f'{intronerate_sequence_directory}/{gene_name}_introns.fasta')
        shutil.move(f'{intronerate_processing_directory}/{gene_name}_introns.fasta',
                    f'{intronerate_sequence_directory}')

    # Move the supercontig sequence to the intronerate_sequence_directory:
    if os.path.exists(f'{intronerate_sequence_directory}/{gene_name}_supercontig.fasta'):
        os.remove(f'{intronerate_sequence_directory}/{gene_name}_supercontig.fasta')
    shutil.move(f'{intronerate_processing_directory}/{gene_name}_supercontig.fasta',
                f'{intronerate_sequence_directory}')


def parse_exonerate_and_get_stitched_contig(exonerate_text_output, query_file, paralog_warning_min_length_percentage,
                                            thresh, logger, prefix, discordant_cutoff, edit_distance, bbmap_subfilter,
                                            bbmap_memory, bbmap_threads, interleaved_fasta_file, no_stitched_contig,
                                            spades_assembly_dict, depth_multiplier, keep_intermediate_files):
    """
    => Parses the C4 alignment text output of Exonerate using BioPython SearchIO.
    => Generates paralog warning and fasta files.
    => Generates stitched contig (or single hit) fasta files for nucleotide and amino-acid sequences.
    => Performs a stitched contig chimera test if file of R1/R2 interleaved reads is present and no_stitched_contig is
    False.

    :param str exonerate_text_output: path to the results text file output by Exonerate (--showalignment yes)
    :param str query_file: path to the protein query fasta file.
    :param float paralog_warning_min_length_percentage: percentage coverage of query required for paralog warning.
    :param int thresh: minimum percentage similarity threshold used to filter Exonerate hits.
    :param logger logger: a logger object
    :param str prefix: path to gene/sample folder e.g. gene001/sampleID
    :param int discordant_cutoff: number of discordant read pairs for a stitched contig to be flagged as chimeric
    :param int edit_distance: edit distance threshold for identifying discordant read pairs
    :param int bbmap_subfilter: ban bbmap.sh alignments with more than this many substitutions
    :param int bbmap_memory: GB of RAM to use for bbmap.sh
    :param int bbmap_threads: number of threads to use for bbmap.sh
    :param None, str interleaved_fasta_file: path to the file of interleaved R1 and R2 fasta seqs, if present
    :param bool no_stitched_contig: if True, return the longest Exonerate hit only
    :param dict spades_assembly_dict: a dictionary of raw SPAdes contigs
    :param int depth_multiplier: assign long paralog as main if coverage depth <depth_multiplier> time other paralogs
    :param bool keep_intermediate_files: if True, keep intermediate files from stitched contig
    :return __main__.Exonerate: instance of the class Exonerate for a given gene
    """

    # Note that the BioPython SearchIO module still has some bugs when Exonerate returns an alignment with directly
    #  abutting introns, or an alignment that ends on an intron. Will submit bug-fixes, but in the meantime catch any
    #  errors and move on:

    try:
        exonerate_hits_from_alignment = list(SearchIO.parse(exonerate_text_output, 'exonerate-text'))
    except AttributeError as e:
        logger.debug(f'Error parsing Exonerate results file for prefix {prefix}\nError is: {e} ')
        return None

    logger.debug(f'no_stitched_contig is: {no_stitched_contig}')

    exonerate_result = Exonerate(searchio_object=exonerate_hits_from_alignment,
                                 query_file=query_file,
                                 paralog_warning_min_length_percentage=paralog_warning_min_length_percentage,
                                 thresh=thresh,
                                 logger=logger,
                                 prefix=prefix,
                                 discordant_cutoff=discordant_cutoff,
                                 edit_distance=edit_distance,
                                 bbmap_subfilter=bbmap_subfilter,
                                 bbmap_memory=bbmap_memory,
                                 bbmap_threads=bbmap_threads,
                                 interleaved_fasta_file=interleaved_fasta_file,
                                 no_stitched_contig=no_stitched_contig,
                                 spades_assembly_dict=spades_assembly_dict,
                                 depth_multiplier=depth_multiplier,
                                 keep_intermediate_files=keep_intermediate_files)

    logger.debug(exonerate_result)

    if not exonerate_result.hits_filtered_by_pct_similarity_dict:  # i.e. no hits left after filtering via pct ID
        return None

    # if exonerate_result.long_paralogs_dict:  # i.e. there are long paralogs recovered
    exonerate_result.write_long_paralogs_and_warnings_to_file()

    if no_stitched_contig:
        exonerate_result.write_no_stitched_contig()
    else:
        exonerate_result.write_stitched_contig_to_file()

    if keep_intermediate_files:
        exonerate_result.write_trimmed_stitched_contig_hits_to_file()
    exonerate_result.write_exonerate_stats_file()

    return exonerate_result


class Exonerate(object):
    """
    Class to parse Exonerate results (SearchIO object) for a given gene. Returns an Exonerate object.
    """

    def __init__(self,
                 searchio_object,
                 query_file=None,
                 paralog_warning_min_length_percentage=0.75,
                 thresh=55,
                 logger=None,
                 prefix=None,
                 discordant_cutoff=5,
                 edit_distance=5,
                 bbmap_subfilter=7,
                 bbmap_memory=1,
                 bbmap_threads=1,
                 interleaved_fasta_file=None,
                 no_stitched_contig=False,
                 spades_assembly_dict=None,
                 depth_multiplier=10,
                 keep_intermediate_files=False):
        """
        Initialises class attributes.

        :param list searchio_object: list returned by parsing Exonerate output (--showalignment yes) with SearchIO.parse
        :param str query_file: path to the protein query fasta file
        :param float paralog_warning_min_length_percentage: percentage coverage of query required for paralog warning
        :param int thresh: minimum percentage similarity threshold used to filter Exonerate hits
        :param logging.RootLogger logger: a logger object
        :param str prefix: path to gene/sample folder e.g. gene001/sampleID
        :param int discordant_cutoff: number of discordant read pairs for a stitched contig to be flagged as chimeric
        :param int edit_distance: edit distance threshold for identifying discordant read pairs
        :param int bbmap_subfilter: ban bbmap.sh alignments with more than this many substitutions
        :param int bbmap_memory: GB of RAM to use for bbmap.sh
        :param int bbmap_threads: number of threads to use for bbmap.sh
        :param str interleaved_fasta_file: path to the file of interleaved R1 and R2 fasta seqs, if present
        :param bool no_stitched_contig: if True, return the longest Exonerate hit only
        :param dict spades_assembly_dict: a dictionary of raw SPAdes contigs
        :param int depth_multiplier: assign long paralog as main if coverage depth <depth_multiplier> other paralogs
        :param bool keep_intermediate_files: if True, keep intermediate files from stitched contig
        """

        if len(searchio_object) != 1:  # This should always be 1 for a single Exonerate query
            raise ValueError(f'searchio_object list is greater than 1!')

        self.exonerate_searchio_alignment = searchio_object
        self.query_id = searchio_object[0].id
        self.query_length = len(SeqIO.read(query_file, 'fasta'))
        self.similarity_threshold = thresh
        self.paralog_warning_by_contig_length_pct = paralog_warning_min_length_percentage
        self.logger = logger
        self.prefix = prefix
        self.no_stitched_contig = no_stitched_contig
        self.interleaved_fasta_file = interleaved_fasta_file
        self.chimera_discordant_cutoff = discordant_cutoff
        self.chimera_edit_distance = edit_distance
        self.chimera_bbmap_subfilter = bbmap_subfilter
        self.chimera_bbmap_memory = bbmap_memory
        self.chimera_bbmap_threads = bbmap_threads
        self.spades_assembly_dict = spades_assembly_dict
        self.depth_multiplier = depth_multiplier
        self.keep_intermediate_files = keep_intermediate_files
        self.hits_filtered_by_pct_similarity_dict = self._parse_searchio_object()
        self.hits_subsumed_hits_removed_dict = self._remove_subsumed_hits()
        self.hits_subsumed_hits_removed_overlaps_trimmed_dict = self._trim_overlapping_hits()
        self.long_paralogs_dict = self._recover_long_paralogs()
        self.paralog_warning_by_contig_depth = self._paralog_warning_by_contig_depth()
        self.stitched_contig_seqrecord = self._create_stitched_contig()
        self.stitched_contig_hit_ranges = self._get_stitched_contig_hit_ranges()
        # only generate a no_stitched_contig_seqrecord (and write report) if no_stitched_contig is True:
        if self.no_stitched_contig:
            self.no_stitched_contig_seqrecord = self._no_stitched_contig()
        else:
            self.no_stitched_contig_seqrecord = None

        # Only perform chimera test if stitched contigs are being created AND interleaved_fasta_file is not None AND a
        # multi-hit stitched contig has been created:
        if self.hits_filtered_by_pct_similarity_dict and not self.no_stitched_contig and interleaved_fasta_file and \
                not self.stitched_contig_seqrecord.description == 'single_hit':
            self.chimera_warning_bool = self._stitched_contig_chimera_warning()
        else:
            self.chimera_warning_bool = None

    def _parse_searchio_object(self):
        """
        Parses the object returned by BioPython SearchIO.parse.
        => Calculates query-vs-hit similarity scores for each hit, and filters hit based on a given threshold
        => Sorts similarity-filtered hsps by start position in the protein query
        => Populates a dict of dicts for each hit, with hitname: {key:value hit data}

        :return collections.defaultdict filtered_by_similarity_hsps_dict: dict of dicts for each hit
        """

        single_exonerate_qresult = self.exonerate_searchio_alignment[0]
        filtered_hsps = []

        # Calculate hsp similarity and filter hsps via a given threshold similarity percentage:
        for hsp in single_exonerate_qresult.hsps:
            similarity_count_total = 0
            similarity_count = 0
            for alignment in hsp.aln_annotation_all:
                for triplet in alignment['similarity']:
                    if triplet == '|||':
                        similarity_count += 1
                    similarity_count_total += 1
            similarity = f'{similarity_count / similarity_count_total * 100:.2f}'
            if float(similarity) > self.similarity_threshold:
                filtered_hsps.append([hsp, float(similarity)])  # add list with hsp object and similarity percentage

        self.logger.debug(f'Number of HSPs before filtering via percent_similarity:'
                          f' {len(single_exonerate_qresult.hsps)}')
        self.logger.debug(f'Number of HSPs after filtering via percent_similarity: {len(filtered_hsps)}')

        if len(filtered_hsps) == 0:
            self.logger.debug(f'Number of HSPs after filtering via percent_similarity is zero')
            return None

        filtered_by_similarity_hsps_dict = defaultdict(dict)  # dict of dicts for each filtered hsp

        for filtered_hsp in sorted(filtered_hsps, key=lambda x: x[0].query_start):  # sort hsps by query start location
            spades_contig_depth = float(filtered_hsp[0].hit_id.split('_')[-1])  # dependant on SPAdes header output!
            query_range = filtered_hsp[0].query_range
            query_range_all = filtered_hsp[0].query_range_all
            hit_range = filtered_hsp[0].hit_range
            hit_range_all = filtered_hsp[0].hit_range_all
            hit_inter_ranges = filtered_hsp[0].hit_inter_ranges
            hit_similarity = filtered_hsp[1]
            hsp_hit_strand_all = filtered_hsp[0].hit_strand_all
            assert len(set(hsp_hit_strand_all)) == 1  # Check that all HSP fragments are on the same strand
            hsp_hit_strand = next(iter(set(hsp_hit_strand_all)))

            # Set a unique hit name for cases where there's >1 hit for a single SPAdes contig:
            unique_hit_name = f'{filtered_hsp[0].hit_id},{self.query_id},{query_range[0]},{query_range[1]}' \
                              f',{hit_similarity},({hsp_hit_strand}),{hit_range[0]},{hit_range[1]}'

            # Get hit hsp nucleotide sequence (exonic regions only):
            concatenated_hsp_alignment_seqs = []
            for alignment in filtered_hsp[0].aln_annotation_all:
                alignment_seq = ''.join(alignment['hit_annotation'])
                concatenated_hsp_alignment_seqs.append(alignment_seq)
            hit_seq = SeqRecord(id=unique_hit_name, name=unique_hit_name, description=unique_hit_name,
                                seq=Seq(''.join(concatenated_hsp_alignment_seqs)).ungap(gap='-'))

            # Populate nested dictionary for hsp:
            filtered_by_similarity_hsps_dict[unique_hit_name]['query_range'] = query_range
            filtered_by_similarity_hsps_dict[unique_hit_name]['query_range_all'] = query_range_all
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_range'] = hit_range
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_range_all'] = hit_range_all
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_inter_ranges'] = hit_inter_ranges
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_strand'] = hsp_hit_strand
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_sequence'] = hit_seq
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_spades_contig_depth'] = spades_contig_depth
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_similarity'] = hit_similarity

        return filtered_by_similarity_hsps_dict

    def _recover_long_paralogs(self):
        """
        Determines whether there are multiple Exonerate hits that are >75% the length of the query protein. If so,
        the 'main' sequence is selected based first on SPAdes contig depth, then (if not definitive) based on
        similarity to the query sequence.

        :return None, collections.defaultdict: paralog_dicts_by_depth OR paralog_dicts_by_percent_similarity
        """

        if not self.hits_filtered_by_pct_similarity_dict or len(self.hits_filtered_by_pct_similarity_dict) == 1:  # i.e.
            # only one hit, so no paralogs
            return None

        paralog_dicts = defaultdict(dict)
        for hit_dict_key, hit_dict_values in self.hits_filtered_by_pct_similarity_dict.items():
            total_hit_vs_query_coverage_length = 0
            for query_range in hit_dict_values['query_range_all']:
                assert query_range[1] >= query_range[0]  # >= required as frameshift can result in e.g. (95, 95)
                hit_vs_query_coverage = query_range[1] - query_range[0]
                total_hit_vs_query_coverage_length += hit_vs_query_coverage

            # Check if hit is longer than threshold:
            if total_hit_vs_query_coverage_length / self.query_length > \
                    self.paralog_warning_by_contig_length_pct:
                paralog_dicts[hit_dict_key] = copy.deepcopy(hit_dict_values)

        if len(paralog_dicts) <= 1:  # i.e. multiple long paralogs not found
            return None

        # Assign '*.main' seq via depth or else percent similarity (written to SeqRecord description):
        paralog_dicts_by_depth = self._best_long_paralog_by_depth(paralog_dicts)
        if paralog_dicts_by_depth:
            return paralog_dicts_by_depth
        paralog_dicts_by_percent_similarity = self._best_long_paralog_by_percent_similarity(paralog_dicts)
        if paralog_dicts_by_percent_similarity:
            return paralog_dicts_by_percent_similarity
        else:
            raise ValueError(f'Issue naming paralogs, please check!')

    def _best_long_paralog_by_depth(self, paralog_dicts):
        """
        Checks if one of the paralogs in the dict paralog_dicts has a SPAdes coverage value >=10 times (default) all
        other paralogs. If so, annotates the high-depth paralog SeqRecord description. If not, returns None.

        :param collections.defaultdict paralog_dicts: dictionary of dictionaries; hitname: {key:value hit data}
        :return NoneType, collections.defaultdict: paralog_dicts with seq.description 'paralog_main_by_depth'
        """

        all_names_and_depths = []
        for paralog_name, paralog_data_dict in paralog_dicts.items():
            all_names_and_depths.append((paralog_name, paralog_data_dict['hit_spades_contig_depth']))
        all_names_and_depths.sort(reverse=True, key=itemgetter(1))
        max_depth_paralog_name = all_names_and_depths[0][0]
        max_depth = all_names_and_depths[0][1]  # single float
        remaining_depths = [item[1] for item in all_names_and_depths[1:]]  # list of one or more floats
        if len(remaining_depths) == 0:  # i.e. there's only one paralog sequence
            paralog_dicts[max_depth_paralog_name]['hit_sequence'].description = 'paralog_main_by_default'
            return paralog_dicts

        depth_threshold = max_depth / self.depth_multiplier  # default depth_multiplier is 10
        try:
            assert all(depth <= depth_threshold for depth in remaining_depths)
            paralog_dicts[max_depth_paralog_name]['hit_sequence'].description = 'paralog_main_by_depth'
            return paralog_dicts
        except AssertionError:
            return None

    @staticmethod
    def _best_long_paralog_by_percent_similarity(paralog_dicts):
        """
        Checks which of the paralogs in the dict paralog_dicts has the highest percentage similarity to the query
        sequence, and annotates the high-similarity paralog SeqRecord description. If there are multiple paralogs
        with an equal high similarity percentage, the annotated paralog will be the first in the list return by Python's
        list.sort() method.

        :param collections.defaultdict paralog_dicts: dictionary of dictionaries; hitname: {key:value hit data}
        :return collections.defaultdict: paralog_dicts, with highest similarity paralog annotated in seq.description
        """

        all_names_and_percent_ids = []
        for paralog_name, paralog_data_dict in paralog_dicts.items():
            all_names_and_percent_ids.append((paralog_name, paralog_data_dict['hit_similarity']))
        all_names_and_percent_ids.sort(reverse=True, key=itemgetter(1))
        max_percent_similarity_paralog_name = all_names_and_percent_ids[0][0]
        paralog_dicts[max_percent_similarity_paralog_name]['hit_sequence'].description = 'paralog_main_by_percent_id'
        return paralog_dicts

    def write_long_paralogs_and_warnings_to_file(self):
        """
        => Renames long paralog sequences for writing to fasta file, using the suffix'.main' for the 'main' selected
        paralog, and then incrementing suffixes for the remaining paralogs ('*.0', '*.1', ...).
        => Writes a fasta file of long paralog sequences.
        => Writes a *.txt file with paralog warnings for long paralogs.

        :return NoneType: no explicit return
        """

        # Create a folder to write long paralog fasta files:
        sample_name = os.path.split(self.prefix)[-1]
        gene_name = os.path.split(self.prefix)[-2]
        paralogs_folder = f'{self.prefix}/paralogs'

        paralog_count = 0
        paralog_seqs_to_write = []
        paralog_warning_strings = []

        if self.long_paralogs_dict:
            if not os.path.exists(paralogs_folder):  # only make directory if long paralogs are present
                os.mkdir(paralogs_folder)

            for paralog_name, paralog_data_dict in self.long_paralogs_dict.items():
                paralog_warning_strings.append(f'{self.query_id}\t{paralog_name}')

                # Rename paralog sequences for writing to fasta file:
                if re.search('main', paralog_data_dict['hit_sequence'].description):
                    paralog_name_to_write = f'{sample_name}.main'
                    paralog_data_dict['hit_sequence'].name = paralog_name_to_write
                    paralog_data_dict['hit_sequence'].id = paralog_name_to_write
                    paralog_data_dict['hit_sequence'].description = paralog_name  # i.e. original SPAdes name
                    paralog_seqs_to_write.append(paralog_data_dict['hit_sequence'])
                else:
                    paralog_name_to_write = f'{sample_name}.{paralog_count}'
                    paralog_data_dict['hit_sequence'].name = paralog_name_to_write
                    paralog_data_dict['hit_sequence'].id = paralog_name_to_write
                    paralog_data_dict['hit_sequence'].description = paralog_name
                    paralog_seqs_to_write.append(paralog_data_dict['hit_sequence'])
                    paralog_count += 1

            # Write fasta file of long paralog sequences:
            with open(f'{paralogs_folder}/{gene_name}_paralogs.fasta', 'w') as paralogs_fasta_handle:
                SeqIO.write(paralog_seqs_to_write, paralogs_fasta_handle, 'fasta')

            # Write *.txt file with paralog warnings for long paralogs
            with open(f'{self.prefix}/paralog_warning_long.txt', 'w') as paralogs_warning_long_handle:
                for warning_string in paralog_warning_strings:
                    paralogs_warning_long_handle.write(f'{warning_string}\n')

        # Write *.txt file with paralog warnings based on contig depth across query:
        with open(f'{self.prefix}/paralog_warning_by_contig_depth.txt', 'w') as paralogs_warning_short_handle:
            paralogs_warning_short_handle.write(f'Sample {sample_name}, gene {gene_name}, paralog warning via contig '
                                                f'depth across query protein,'
                                                f' {self._paralog_warning_by_contig_depth()}\n')

    def _paralog_warning_by_contig_depth(self):
        """
        Takes a dictionary of hits filtered by percent ID, and calculates contig coverage across the length of the
        query protein. If coverage is >1 for a given percentage threshold of the query length, return True.

        :return bool: True if hit coverage is >1 for a given percentage length of the query
        """

        if not self.hits_filtered_by_pct_similarity_dict:
            return None

        # Get tuples of the query ranges for all similarity-filtered hits:
        hit_vs_query_ranges_all = []
        for hit, hit_dict_values in self.hits_filtered_by_pct_similarity_dict.items():
            for query_range in hit_dict_values['query_range_all']:
                assert query_range[1] >= query_range[0]  # >= required as frameshift can result in e.g. (95, 95)
                hit_vs_query_range = (query_range[0], query_range[1])
                hit_vs_query_ranges_all.append(hit_vs_query_range)

        # Get a count of hit coverage for each position in the query sequence:
        query_coverage = []
        for position in range(0, self.query_length):
            coverage = 0
            for hit_vs_query_range in hit_vs_query_ranges_all:
                if position in range(hit_vs_query_range[0], hit_vs_query_range[1]):
                    coverage += 1
            query_coverage.append(coverage)

        # If the hit coverage is >1 for a given percentage length of the query, return True:
        pct_query_cov_greater_than_one = len([depth for depth in query_coverage if depth > 1]) / len(query_coverage)
        if pct_query_cov_greater_than_one > self.paralog_warning_by_contig_length_pct:
            return True
        else:
            return False

    def _remove_subsumed_hits(self):
        """
        Takes a dictionary of hits that have been filtered via similarity to the query, and removes any hit that has
        a query range that completely subsumes (i.e. encompasses) another hit. If two hits have an identical query
        range (and are not themselves subsumed by another longer hit) the hit with the highest similarity is retained.

        :return collections.defaultdict: exonerate_hits_filtered_no_subsumed
        """

        if not self.hits_filtered_by_pct_similarity_dict:
            return None

        exonerate_hits_filtered_no_subsumed = copy.deepcopy(self.hits_filtered_by_pct_similarity_dict)
        hit_comparisons = itertools.permutations(self.hits_filtered_by_pct_similarity_dict, 2)

        seqs_removed = []  # Use to avoid trying combinations with previously removed hits
        hits_with_identical_range_and_similarity_dict = defaultdict(set)
        for hit_pair in hit_comparisons:
            to_remove = None

            hit_1_query_range = self.hits_filtered_by_pct_similarity_dict[hit_pair[0]]['query_range']
            hit_1_similarity = self.hits_filtered_by_pct_similarity_dict[hit_pair[0]]['hit_similarity']
            hit_2_query_range = self.hits_filtered_by_pct_similarity_dict[hit_pair[1]]['query_range']
            hit_2_similarity = self.hits_filtered_by_pct_similarity_dict[hit_pair[1]]['hit_similarity']

            if hit_1_query_range[0] < hit_2_query_range[0] and hit_1_query_range[1] > hit_2_query_range[1]:
                to_remove = hit_pair[1]
                if to_remove in seqs_removed:
                    continue
                seqs_removed.append(to_remove)
            elif hit_1_query_range[0] < hit_2_query_range[0] and hit_1_query_range[1] == hit_2_query_range[1]:
                to_remove = hit_pair[1]
                if to_remove in seqs_removed:
                    continue
                seqs_removed.append(to_remove)
            elif hit_1_query_range[0] == hit_2_query_range[0] and hit_1_query_range[1] > hit_2_query_range[1]:
                to_remove = hit_pair[1]
                if to_remove in seqs_removed:
                    continue
                seqs_removed.append(to_remove)
            else:
                if hit_1_query_range[0] == hit_2_query_range[0] and hit_1_query_range[1] == hit_2_query_range[1]:
                    if hit_1_similarity > hit_2_similarity:
                        to_remove = hit_pair[1]
                    elif hit_1_similarity < hit_2_similarity:
                        to_remove = hit_pair[0]
                    elif hit_1_similarity == hit_2_similarity:  # capture these in a dict and select one for each range
                        # TODO: compare SPAdes depth to select a hit sequence
                        hits_with_identical_range_and_similarity_dict[hit_1_query_range].add(hit_pair[0])
                        hits_with_identical_range_and_similarity_dict[hit_1_query_range].add(hit_pair[1])
                        continue

                    if to_remove in seqs_removed:
                        self.logger.debug(f'to_remove {to_remove} is already in seqs_removed: {seqs_removed}')
                        continue
                    seqs_removed.append(to_remove)

            # If hit not already removed, remove it from the dict:
            if to_remove:
                try:
                    del exonerate_hits_filtered_no_subsumed[to_remove]
                except KeyError:
                    self.logger.debug(f'hit {to_remove} already removed from dict')

        # Select a single sequence from each range in hits_with_identical_range_and_similarity_dict:
        if len(hits_with_identical_range_and_similarity_dict) != 0:
            self.logger.debug(f'Gene has hits with identical query ranges and similarities; selecting one hit for '
                              f'each range')
            self.logger.debug(f'Dictionary hits_with_identical_range_and_similarity_dict is:'
                              f' {hits_with_identical_range_and_similarity_dict}')
            for query_range, hits in hits_with_identical_range_and_similarity_dict.items():
                to_remove = list(hits)[0]  # arbitrarily remove first hit if range and similarity are the same
                try:
                    del exonerate_hits_filtered_no_subsumed[to_remove]
                except KeyError:
                    self.logger.debug(f'hit {to_remove} already removed from dict')

        return exonerate_hits_filtered_no_subsumed

    def _trim_overlapping_hits(self):  # for constructing the coding-seq-only stitched contig via
        # _create_stitched_contig()
        """
        => Takes a dictionary of hits that has been filtered via hit similarity to the query, and has had subsumed hits
        removed. If any of the remaining hits have overlaps in query ranges, the 3' end of the left hit is trimmed to
        remove overlap sequence.
        => If there are any gaps between hit pairs with respect to the query sequence, pads the 3' end of the left hit
        with a corresponding number of ends Ns (query range gaps converted from #amino-acids to #nucleotides).
        => For any trimmed hit, annotate the description of the corresponding SeqRecord with '3prime overlap trimmed by:
        <int>'.

        :return collections.defaultdict: exonerate_hits_subsumed_and_trimmed_dict
        """

        if not self.hits_filtered_by_pct_similarity_dict:
            return None

        if len(self.hits_subsumed_hits_removed_dict) == 1:  # i.e. single hit remaining from previous filtering
            for key, value in self.hits_subsumed_hits_removed_dict.items():
                value['hit_sequence'].description = f'Single hit after filtering: N/A'
            return self.hits_subsumed_hits_removed_dict

        # Don't overwrite self.hits_subsumed_hits_removed_dict:
        exonerate_hits_subsumed_hits_removed_copy = copy.deepcopy(self.hits_subsumed_hits_removed_dict)

        exonerate_hits_subsumed_and_trimmed_dict = {}
        for pair in list(pairwise_longest(exonerate_hits_subsumed_hits_removed_copy.values())):
            if pair[1] is not None:  # as pairwise_longest() will pad a hit without a pair with 'None'
                left_seq_hit_name, right_seq_hit_name = pair[0]['hit_sequence'].id, pair[1]['hit_sequence'].id
                left_seq_query_range, right_seq_query_range = pair[0]['query_range'], pair[1]['query_range']

                # If overlapping hits, always trim the 3' end of the left hit:
                if left_seq_query_range[1] > right_seq_query_range[0]:
                    num_bases_overlap = (left_seq_query_range[1] - right_seq_query_range[0]) * 3
                    bases_to_recover = len(pair[0]['hit_sequence'].seq) - num_bases_overlap
                    pair[0]['hit_sequence'].seq = pair[0]['hit_sequence'].seq[:bases_to_recover]
                    pair[0]['hit_sequence'].description = f'3prime overlap trimmed by: {num_bases_overlap}'
                    exonerate_hits_subsumed_and_trimmed_dict[left_seq_hit_name] = pair[0]

                else:  # if no overlap, add left hit unmodified (3' padded with Ns if gap before next hit):
                    bases_in_gap_between_hits = (right_seq_query_range[0] - left_seq_query_range[1]) * 3
                    pair[0]['hit_sequence'].seq = Seq(f"{pair[0]['hit_sequence'].seq}"
                                                      f"{'N' * bases_in_gap_between_hits}")
                    pair[0]['hit_sequence'].description = f'No overlap: N/A'
                    exonerate_hits_subsumed_and_trimmed_dict[left_seq_hit_name] = pair[0]

            else:  # process final unpaired left hit
                left_seq_hit_name = pair[0]['hit_sequence'].id
                pair[0]['hit_sequence'].description = f'No overlap: N/A'
                exonerate_hits_subsumed_and_trimmed_dict[left_seq_hit_name] = pair[0]

        return exonerate_hits_subsumed_and_trimmed_dict

    def _create_stitched_contig(self):  # Here we're not dealing with intron sequence at all
        """
        Takes a dictionary of filtered, sorted, and trimmed hits, and creates a stitched contig SeqRecord by
        concatenating the corresponding sequences. If only one hit is present, return a SeqRecord of this hit.

        :return Bio.SeqRecord.SeqRecord: no_stitched_contig or stitched_contig, depending on number of hits
        """

        if not self.hits_filtered_by_pct_similarity_dict:
            return None

        sample_name = os.path.split(self.prefix)[-1]
        gene_name = os.path.split(self.prefix)[-2]

        if len(self.hits_subsumed_hits_removed_overlaps_trimmed_dict) == 1:  # i.e. only one hit
            for hit, hit_dict_values in self.hits_subsumed_hits_removed_overlaps_trimmed_dict.items():
                no_stitched_contig_seqrecord = SeqRecord(
                    seq=hit_dict_values['hit_sequence'].seq, id=sample_name, name=sample_name,
                    description='single_hit')

                # Write report file:
                if not self.no_stitched_contig:
                    log_entry = f'{sample_name},{gene_name}, No stitched contig produced. Gene sequence contains a ' \
                                f'single Exonerate hit.'
                    self._write_genes_with_stitched_contig(log_entry)

                return no_stitched_contig_seqrecord

        # If multiple hits:
        stitched_contig_hits = []
        for hit, hit_dict_values in self.hits_subsumed_hits_removed_overlaps_trimmed_dict.items():
            stitched_contig_hits.append(str(hit_dict_values['hit_sequence'].seq))
        num_hits_in_stitched_contig = len(stitched_contig_hits)
        stitched_contig_seqrecord = SeqRecord(seq=Seq(''.join(stitched_contig_hits)), id=sample_name, name=sample_name,
                                              description=f'multi_hit_stitched_contig_comprising_'
                                                          f'{num_hits_in_stitched_contig}_hits')

        # Write report file:
        if not self.no_stitched_contig:
            log_entry = f'{sample_name},{gene_name}, Stitched contig produced. Gene sequence contains more than one ' \
                        f'Exonerate hit.'
            self._write_genes_with_stitched_contig(log_entry)

        return stitched_contig_seqrecord

    @staticmethod
    def convert_coords_revcomp(list_of_range_tuples, raw_spades_contig_length):
        """
        This function takes a list of Exonerate SearchIO hit range tuples for a hit on the negative strand,
        and converts them so that they are consistent with those from a hit on the positive strand. The
        reverse-complemented SPAdes contig can then be processed with the approach used for positive strand
        hits.

        => e.g. for SPAdes contig with length 873 and hit range list [(284, 377), (2, 119)]:
        [(284, 377), (2, 119)] -> [(496, 589), (754, 871)]

        :param list list_of_range_tuples: list of Exonerate SearchIO hit range tuples
        :param int raw_spades_contig_length: integer corresponding to the length of the SPAdes contig
        :return list converted_list: list of hit range tuples converted to positive strand coordinates
        """

        range_tuples_reversed = [tuple(reversed(revcomp_hit_range)) for revcomp_hit_range in
                                 list_of_range_tuples]
        converted_list = [(raw_spades_contig_length - revcomp_hit_range[0], raw_spades_contig_length -
                           revcomp_hit_range[1]) for revcomp_hit_range in range_tuples_reversed]
        return converted_list

    def _get_stitched_contig_hit_ranges(self):
        """
        Returns a dictionary of hit_id:exon coordinates for the full exon-only stitched contig fasta sequence (i.e.
        starting at nucleotide zero and with introns removed).

        :return dict hit_ranges_dict: a dictionary of hit_id: exon coordinates
        """

        if not self.hits_filtered_by_pct_similarity_dict:
            return None

        hit_ranges_dict = defaultdict(list)

        cumulative_hit_length = 0  # track to adjust coordinates of hits to match exon-only stitched contig sequence
        for hit, hit_dict_values in self.hits_subsumed_hits_removed_overlaps_trimmed_dict.items():
            self.logger.debug(f'hit is: {hit}')
            spades_name = hit.split(',')[0]
            raw_spades_contig_length = len(self.spades_assembly_dict[spades_name])
            hit_exonerate_sequence_length = len(hit_dict_values['hit_sequence'].seq)
            hit_range_all = hit_dict_values['hit_range_all']
            hit_inter_ranges = hit_dict_values['hit_inter_ranges']

            if len(hit_inter_ranges) == 0:  # i.e. no introns
                self.logger.debug(f'len(hit_inter_ranges) for hit {hit} is 0, no intron coordinates for chimera test')
                hit_ranges_dict[hit].append('no introns')
                cumulative_hit_length += hit_exonerate_sequence_length
            else:
                if hit_dict_values['hit_strand'] == -1:  # Convert ranges so that they apply to the revcomp contig
                    self.logger.debug(f'hit_inter_ranges before conversion is: {hit_inter_ranges}')
                    hit_inter_ranges = self.convert_coords_revcomp(hit_inter_ranges, raw_spades_contig_length)
                    self.logger.debug(f'hit_inter_ranges after conversion is: {hit_inter_ranges}')
                    self.logger.debug(f'hit_range_all before conversion is: {hit_range_all}')
                    hit_range_all = self.convert_coords_revcomp(hit_range_all, raw_spades_contig_length)
                    self.logger.debug(f'hit_range_all after conversion is: {hit_range_all}')
                else:
                    self.logger.debug(f'hit_inter_ranges is: {hit_inter_ranges}')
                    self.logger.debug(f'hit_range_all is: {hit_range_all}')

                # Adjust range coordinates so that they start at zero (i.e. first nucleotide position of the Exonerate
                # fasta sequence for this hit:
                hit_start_coordinate = hit_range_all[0][0]
                hit_range_all_start_at_zero = [(item[0] - hit_start_coordinate, item[1] - hit_start_coordinate)
                                               for item in hit_range_all]
                self.logger.debug(f'hit_range_all_start_at_zero is: {hit_range_all_start_at_zero}')

                # Adjust range coordinates to remove intron lengths:
                cumulative_intron_length = 0
                hit_ranges_dict[hit].append(hit_range_all_start_at_zero[0])  # add the first exon
                for pair in list(pairwise_longest(hit_range_all_start_at_zero)):
                    if pair[1] is not None:  # as pairwise_longest() will pad a range tuple without a pair with 'None'
                        intron_length = pair[1][0] - pair[0][1]
                        cumulative_intron_length += intron_length
                        no_intron_coordinates = (pair[1][0] - cumulative_intron_length,
                                                 pair[1][1] - cumulative_intron_length)
                        no_intron_coordinates = tuple(no_intron_coordinates)
                        hit_ranges_dict[hit].append(no_intron_coordinates)

                # Adjust hit ranges to account for previous contigs in the stitched_contig:
                hit_ranges_adjusted = [(int(item[0]) + cumulative_hit_length, int(item[1]) + cumulative_hit_length) for
                                       item in hit_ranges_dict[hit]]
                hit_ranges_dict[hit] = hit_ranges_adjusted  # replace with adjusted values

        self.logger.debug(f'hit_ranges_dict is: {hit_ranges_dict}')

        return hit_ranges_dict

    def write_trimmed_stitched_contig_hits_to_file(self):
        """
        Writes DNA (suffix '.FNA') and amino-acid (suffix '.FAA') sequences for the filtered, sorted and trimmed
        Exonerate hits to fasta file. Used for debugging.

        :return NoneType: no explicit return
        """

        # Write DNA seqs
        with open(f'{self.prefix}/exonerate_hits_trimmed.FNA', 'w') as fasta_handle_nucl:
            for key, value in self.hits_subsumed_hits_removed_overlaps_trimmed_dict.items():
                SeqIO.write(value['hit_sequence'], fasta_handle_nucl, 'fasta')

        # Write amino-acid seqs
        with open(f'{self.prefix}/exonerate_hits_trimmed.FAA', 'w') as fasta_handle_amino:
            for key, value in self.hits_subsumed_hits_removed_overlaps_trimmed_dict.items():
                seq_id = value['hit_sequence'].id
                name = value['hit_sequence'].name
                description = value['hit_sequence'].description
                translated_seq = value['hit_sequence'].seq.translate()
                translated_seq_seqrecord = SeqRecord(seq=translated_seq, id=seq_id, name=name, description=description)
                SeqIO.write(translated_seq_seqrecord, fasta_handle_amino, 'fasta')

    def write_exonerate_stats_file(self):
        """
        Write a *.csv file with stats from the initial Exonerate run following filtering by similarity,
        and stats for subsequent filtering/trimming steps.

        :return NoneType: no explicit return
        """

        if not self.hits_filtered_by_pct_similarity_dict:
            self.logger.info(f'There are no Exonerate hits remaining after filtering with a '
                             f'{self.similarity_threshold} percent similarity threshold!')
            return

        sample_name = os.path.split(self.prefix)[-1]
        gene_name = os.path.split(self.prefix)[-2]

        headers = ['query_id', 'query_length', 'hit_id', 'query_range_limits', 'query_hsps_ranges',
                   'hit_percent_similarity', 'hit_strand', 'hit_range_limits', 'hit_hsps_ranges',
                   '3-prime_bases_trimmed']

        def nested_dict_iterator(nested_dict):
            """
            Nested function to yield the key and value of a given nested dict.

            :param dict, collections.defaultdict nested_dict:
            :return: yields the key and value of the nested dict
            """

            for dict_key, dict_value in nested_dict.items():
                if isinstance(dict_value, dict):
                    yield dict_key, dict_value

        with open(f'{self.prefix}/exonerate_stats.tsv', 'w') as exonerate_stats_handle:
            headers = '\t'.join(headers)
            exonerate_stats_handle.write(f"\t{headers}\n")

            # Write stats for hits filtered by similarity:
            nested_dict_iter = nested_dict_iterator(self.hits_filtered_by_pct_similarity_dict)
            first_hit_dict_key, first_hit_dict_value = next(nested_dict_iter)  # get first line
            exonerate_stats_handle.write((f"Hits filtered > {self.similarity_threshold} percent similarity"
                                          f"\t{self.query_id}"
                                          f"\t{self.query_length}"
                                          f"\t{str(first_hit_dict_key.split(',')[0])}"
                                          f"\t{first_hit_dict_value['query_range']}"
                                          f"\t{first_hit_dict_value['query_range_all']}"
                                          f"\t{first_hit_dict_value['hit_similarity']}"
                                          f"\t{first_hit_dict_value['hit_strand']}"
                                          f"\t{first_hit_dict_value['hit_range']}"
                                          f"\t{first_hit_dict_value['hit_range_all']}"
                                          f"\tN/A\n"))  # as no trimming performed yet
            for key, value in nested_dict_iter:  # Write remaining lines
                exonerate_stats_handle.write((f"\t{self.query_id}"
                                              f"\t{self.query_length}"
                                              f"\t{str(key.split(',')[0])}"
                                              f"\t{value['query_range']}"
                                              f"\t{value['query_range_all']}"
                                              f"\t{value['hit_similarity']}"
                                              f"\t{value['hit_strand']}"
                                              f"\t{value['hit_range']}"
                                              f"\t{value['hit_range_all']}"
                                              f"\tN/A\n"))

            # Write stats for hits filtered by similarity and with subsumed hits removed:
            nested_dict_iter = nested_dict_iterator(self.hits_subsumed_hits_removed_dict)
            first_hit_dict_key, first_hit_dict_value = next(nested_dict_iter)  # get first line
            exonerate_stats_handle.write((f"Hits with subsumed hits removed"
                                          f"\t{self.query_id}"
                                          f"\t{self.query_length}"
                                          f"\t{str(first_hit_dict_key.split(',')[0])}"
                                          f"\t{first_hit_dict_value['query_range']}"
                                          f"\t{first_hit_dict_value['query_range_all']}"
                                          f"\t{first_hit_dict_value['hit_similarity']}"
                                          f"\t{first_hit_dict_value['hit_strand']}"
                                          f"\t{first_hit_dict_value['hit_range']}"
                                          f"\t{first_hit_dict_value['hit_range_all']}"
                                          f"\tN/A\n"))  # as no trimming performed yet
            for key, value in nested_dict_iter:  # Write remaining lines
                exonerate_stats_handle.write((f"\t{self.query_id}"
                                              f"\t{self.query_length}"
                                              f"\t{str(key.split(',')[0])}"
                                              f"\t{value['query_range']}"
                                              f"\t{value['query_range_all']}"
                                              f"\t{value['hit_similarity']}"
                                              f"\t{value['hit_strand']}"
                                              f"\t{value['hit_range']}"
                                              f"\t{value['hit_range_all']}"
                                              f"\tN/A\n"))

            # Write stats for hits filtered by similarity and subsumed hits removed and trimmed:
            nested_dict_iter = nested_dict_iterator(self.hits_subsumed_hits_removed_overlaps_trimmed_dict)
            first_hit_dict_key, first_hit_dict_value = next(nested_dict_iter)  # get first line
            trimmed_bases = first_hit_dict_value['hit_sequence'].description.split(':')[-1].strip()
            exonerate_stats_handle.write((f"Hits with subsumed hits removed and trimmed"
                                          f"\t{self.query_id}"
                                          f"\t{self.query_length}"
                                          f"\t{str(first_hit_dict_key.split(',')[0])}"
                                          f"\t{first_hit_dict_value['query_range']}"
                                          f"\t{first_hit_dict_value['query_range_all']}"
                                          f"\t{first_hit_dict_value['hit_similarity']}"
                                          f"\t{first_hit_dict_value['hit_strand']}"
                                          f"\t{first_hit_dict_value['hit_range']}"
                                          f"\t{first_hit_dict_value['hit_range_all']}"
                                          f"\t{trimmed_bases}\n"))
            for key, value in nested_dict_iter:  # Write remaining lines
                trimmed_bases = value['hit_sequence'].description.split()[-1]
                exonerate_stats_handle.write((f"\t{self.query_id}"
                                              f"\t{self.query_length}"
                                              f"\t{str(key.split(',')[0])}"
                                              f"\t{value['query_range']}"
                                              f"\t{value['query_range_all']}"
                                              f"\t{value['hit_similarity']}"
                                              f"\t{value['hit_strand']}"
                                              f"\t{value['hit_range']}"
                                              f"\t{value['hit_range_all']}"
                                              f"\t{trimmed_bases}\n"))

            # Write stats for paralogs:
            if self.long_paralogs_dict:
                dict_iter = iter((key, value) for key, value in self.long_paralogs_dict.items())
                first_hit_dict_key, first_hit_dict_value = next(dict_iter)  # get first line
                exonerate_stats_handle.write(
                    (f"Hits corresponding to long paralogs"
                     f"\t{self.query_id}"
                     f"\t{self.query_length}"
                     f"\t{str(first_hit_dict_key.split(',')[0])}"
                     f"\t{first_hit_dict_value['query_range']}"
                     f"\t{first_hit_dict_value['query_range_all']}"
                     f"\t{first_hit_dict_value['hit_similarity']}"
                     f"\t{first_hit_dict_value['hit_strand']}"
                     f"\t{first_hit_dict_value['hit_range']}"
                     f"\t{first_hit_dict_value['hit_range_all']}"
                     f"\tN/A\n"))
                for key, value in dict_iter:  # Write remaining lines
                    exonerate_stats_handle.write((f"\t{self.query_id}"
                                                  f"\t{self.query_length}"
                                                  f"\t{str(key.split(',')[0])}"
                                                  f"\t{value['query_range']}"
                                                  f"\t{value['query_range_all']}"
                                                  f"\t{value['hit_similarity']}"
                                                  f"\t{value['hit_strand']}"
                                                  f"\t{value['hit_range']}"
                                                  f"\t{value['hit_range_all']}"
                                                  f"\tN/A\n"))

    def write_stitched_contig_to_file(self):
        """
        Writes DNA (suffix '.FNA') and amino-acid (suffix '.FAA') sequences for the stitched_contig (or single remaining
        hit sequence) to fasta file.

        :return NoneType: no explicit return
        """

        gene_name = os.path.split(self.prefix)[-2]
        hit_id = self.stitched_contig_seqrecord.id
        name = self.stitched_contig_seqrecord.name
        description = self.stitched_contig_seqrecord.description

        dna_seqrecord_to_write = self.stitched_contig_seqrecord
        amino_acid_seq_to_write = self.stitched_contig_seqrecord.seq.translate()
        amino_acid_seqrecord_to_write = SeqRecord(seq=amino_acid_seq_to_write, id=hit_id, name=name,
                                                  description=description)  # as seq.translate() creates 'unknown'

        with open(f'{self.prefix}/sequences/FNA/{gene_name}.FNA', 'w') as fna_handle:
            SeqIO.write(dna_seqrecord_to_write, fna_handle, 'fasta')

        with open(f'{self.prefix}/sequences/FAA/{gene_name}.FAA', 'w') as faa_handle:
            SeqIO.write(amino_acid_seqrecord_to_write, faa_handle, 'fasta')

    def _stitched_contig_chimera_warning(self):
        """
        Produces a warning (boolean) if there's evidence that the Exonerate hit sequences used to stitch together
        a stitched contig are derived from different paralogs. Can only be performed when R1 and R2 reads are present (
        i.e. it doesn't work with single-end reads).

        => Maps R1 and R2 reads for the sample/gene against the stitched_contig sequence.
        => Counts correctly mapped read pairs where one read maps with 100% length and identity to the stitched contig
        reference, and the other read has a number of substitutions greater than a given threshold. This is
        likely to occur across hit boundaries, where hits are derived from different paralogs.
        => Ignores read pairs unless 1) R1 occurs in a different Exonerate hit contig to R2; 2) they fall entirely
        within exon sequences (i.e. they don't overlap exon-intron boundaries which would cause spurious mismatches)

        :return bool: True is a chimera warning is produced and written to file.
        """

        self.logger.debug(f'self.stitched_contig_hit_ranges is: {self.stitched_contig_hit_ranges}')

        if not self.hits_filtered_by_pct_similarity_dict:
            return None

        sample_name = os.path.split(self.prefix)[-1]
        gene_name = os.path.split(self.prefix)[-2]

        # Write the stitched_contig sequence to fasta file for read mapping via bbmap:
        dna_seqrecord_to_write = self.stitched_contig_seqrecord
        with open(f'{self.prefix}/chimera_test_stitched_contig.fasta', 'w') as chimera_test_handle:
            SeqIO.write(dna_seqrecord_to_write, chimera_test_handle, 'fasta')

        # Map interleaved R1 and R2 reads against the stitched_contig sequence:
        bbmap_command = f'bbmap.sh ' \
                        f'-Xmx{self.chimera_bbmap_memory}g ' \
                        f'-t={self.chimera_bbmap_threads} ' \
                        f'ref={self.prefix}/chimera_test_stitched_contig.fasta ' \
                        f'in={self.interleaved_fasta_file} ' \
                        f'out={self.prefix}/chimera_test_stitched_contig.sam ' \
                        f'interleaved=t ' \
                        f'pairedonly=t ' \
                        f'mappedonly=t ' \
                        f'maxindel=0 ' \
                        f'strictmaxindel=t ' \
                        f'nodisk=t ' \
                        f'subfilter={self.chimera_bbmap_subfilter} ' \
                        f'ambiguous=toss'

        try:
            result = subprocess.run(bbmap_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True, check=True)
            self.logger.debug(f'bbmap_command check_returncode() is: {result.check_returncode()}')
            # self.logger.debug(f'bbmap_command stdout is: {result.stdout}')
            # self.logger.debug(f'bbmap_command stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            self.logger.debug(f'bbmap_command FAILED. Output is: {exc}')
            self.logger.debug(f'bbmap_command stdout is: {exc.stdout}')
            self.logger.debug(f'bbmap_command stderr is: {exc.stderr}')

        # Get a list of individual contig ranges within the stitched_contig (dna_seqrecord_to_write):
        hits_processed = []
        individual_contig_ranges_in_stitched_contig = []
        for hit, hit_data_dict in self.hits_subsumed_hits_removed_overlaps_trimmed_dict.items():
            spades_contig_name = hit.split(',')[0]
            if spades_contig_name in hits_processed:
                # TODO I think this can be removed, but keeping to catch issue so I can test it:
                raise ValueError(f'Chimera test: contig {spades_contig_name} already processed for sample'
                                 f' {sample_name}, gene {gene_name} ')
            contig_length = len(hit_data_dict['hit_sequence'])
            if len(individual_contig_ranges_in_stitched_contig) == 0:  # i.e. it's the first contig
                contig_range_tuple = (hit, 0, contig_length)
                individual_contig_ranges_in_stitched_contig.append(contig_range_tuple)
            else:
                contig_start_coordinate = individual_contig_ranges_in_stitched_contig[-1][-1]  # start after last range
                contig_end_coordinate = contig_start_coordinate + contig_length
                contig_range_tuple = (hit, contig_start_coordinate, contig_end_coordinate)
                individual_contig_ranges_in_stitched_contig.append(contig_range_tuple)

        self.logger.debug(f'individual_contig_ranges_in_stitched_contig is:'
                          f' {individual_contig_ranges_in_stitched_contig}')

        # Check that the max range value corresponds to the length of the stitched contig sequence:
        assert individual_contig_ranges_in_stitched_contig[-1][-1] == len(self.stitched_contig_seqrecord)

        # Parse the sam file produced by bbmap.sh:
        samfile_reads = []
        with open(f'{self.prefix}/chimera_test_stitched_contig.sam') as samfile:
            lines = samfile.readlines()
            for line in lines:
                if not line.startswith('@'):
                    samfile_reads.append(line)

        if not self.keep_intermediate_files:
            os.remove(f'{self.prefix}/chimera_test_stitched_contig.sam')

        # Count the number of discordant read pairs and recover R1 and R2 in a list for further filtering:
        discordant_read_list = []
        for forward, reverse in grouped(samfile_reads, 2):
            forward_edit_distance = (forward.split('\t')[11]).split(':')[2]
            reverse_edit_distance = (reverse.split('\t')[11]).split(':')[2]
            if int(forward_edit_distance) == 0 and int(reverse_edit_distance) >= self.chimera_edit_distance:
                discordant_read_list.append(forward)
                discordant_read_list.append(reverse)
            elif int(reverse_edit_distance) == 0 and int(forward_edit_distance) >= self.chimera_edit_distance:
                discordant_read_list.append(forward)
                discordant_read_list.append(reverse)

        self.logger.debug(f'There are {int(len(discordant_read_list) / 2)} discordant read pairs prior to range '
                          f'filtering')

        # Filter discordant read pairs to remove any that 1) don't have each read mapping to different Exonerate hit
        # contigs; 2) don't fall entirely withing exon sequences (i.e. they overlap exon-intron boundaries which
        # would cause spurious mismatches); 2) don't have each read mapping to different Exonerate hit contigs:
        discordant_read_list_pass_filtering = []
        for forward, reverse in grouped(discordant_read_list, 2):
            forward_start_coordinate = int(forward.split('\t')[3]) - 1  # SAM uses 1-based, adjust to zero-based
            forward_seq_len = len(forward.split('\t')[9])
            forward_end_coordinate = forward_start_coordinate + forward_seq_len
            reverse_start_coordinate = int(reverse.split('\t')[3]) - 1  # SAM uses 1-based, adjust to zero-based
            reverse_seq_len = len(reverse.split('\t')[9])
            reverse_end_coordinate = reverse_start_coordinate + reverse_seq_len

            # Check if each read falls within a different contig:
            forward_enclosing_contig_range = None
            reverse_enclosing_contig_range = None

            for range_tuple in individual_contig_ranges_in_stitched_contig:
                if forward_start_coordinate >= range_tuple[1] and forward_end_coordinate <= range_tuple[2]:
                    forward_enclosing_contig_name = range_tuple[0]
                    forward_enclosing_contig_name_spades_only = range_tuple[0].split(',')[0]
                    forward_enclosing_contig_range = range_tuple
                if reverse_start_coordinate >= range_tuple[1] and reverse_end_coordinate <= range_tuple[2]:
                    reverse_enclosing_contig_name = range_tuple[0]
                    reverse_enclosing_contig_name_spades_only = range_tuple[0].split(',')[0]
                    reverse_enclosing_contig_range = range_tuple
            self.logger.debug(f'forward_enclosing_contig_range is: {forward_enclosing_contig_range}')
            self.logger.debug(f'reverse_enclosing_contig_range is: {reverse_enclosing_contig_range}')

            # Make sure each read is assigned to a contig range, or skip if it overlaps contig join:
            if not forward_enclosing_contig_range or not reverse_enclosing_contig_range:
                self.logger.debug(f'One or both read pairs do not map within the range of a single contig (i.e. they '
                                  f'overlap a coordinate where two contigs have been concatenated). Skipping read pair')
                continue

            forward_contig_exon_ranges = self.stitched_contig_hit_ranges[forward_enclosing_contig_name]
            reverse_contig_exon_ranges = self.stitched_contig_hit_ranges[reverse_enclosing_contig_name]

            self.logger.debug(f'forward_contig_exon_ranges is: {forward_contig_exon_ranges}')
            self.logger.debug(f'reverse_contig_exon_ranges is: {reverse_contig_exon_ranges}')

            # Check if the SPAdes only name each contig are the same:
            if forward_enclosing_contig_name_spades_only == reverse_enclosing_contig_name_spades_only:
                self.logger.debug(f'Both reads occur in the same SPAdes contig. Skipping read pair!')
                continue
            else:  # If contig is not the same, check that reads don't overlap exon-intron boundaries within each contig
                self.logger.debug(f'Reads occur in different contigs. Check for exon-intron overlaps...')
                forward_occurs_in_exons_only = False
                reverse_occurs_in_exons_only = False
                for range_tuple in forward_contig_exon_ranges:
                    if range_tuple == 'no introns':
                        forward_occurs_in_exons_only = True
                        break
                    if forward_start_coordinate >= range_tuple[0] and forward_end_coordinate <= range_tuple[1]:
                        forward_occurs_in_exons_only = True
                        break
                for range_tuple in reverse_contig_exon_ranges:
                    if range_tuple == 'no introns':
                        reverse_occurs_in_exons_only = True
                        break
                    if reverse_start_coordinate >= range_tuple[0] and reverse_end_coordinate <= range_tuple[1]:
                        reverse_occurs_in_exons_only = True
                        break

                # Check both reads occur in exons only, and skip if not:
                if forward_occurs_in_exons_only and reverse_occurs_in_exons_only:
                    discordant_read_list_pass_filtering.append(forward)
                    discordant_read_list_pass_filtering.append(reverse)
                else:
                    self.logger.debug(f'One or both read pairs overlap with exon-intron boundaries - skipping read '
                                      f'pair')

            # If there are discordant read pairs passing filtering, write them to a SAM file, and write to log:
            if discordant_read_list_pass_filtering:
                number_of_discordant_read_pairs_passing_filtering = int(len(discordant_read_list_pass_filtering) / 2)
                with open(f'{self.prefix}/chimera_test_diagnostic_reads.sam', 'w') as diagnostic_reads:
                    for forward_read, reverse_read in grouped(discordant_read_list_pass_filtering, 2):
                        diagnostic_reads.write(forward_read)
                        diagnostic_reads.write(reverse_read)

                if number_of_discordant_read_pairs_passing_filtering > self.chimera_discordant_cutoff:
                    # Write report file for gene
                    with open(f'{self.prefix}/putative_chimeric_stitched_contig.csv', 'w') as \
                            discordant_stitched_contig_reportfile:
                        log_entry = f'{sample_name},{gene_name}, Chimera WARNING for stitched_contig. Sequence may ' \
                                    f'be derived from multiple paralogs.'
                        discordant_stitched_contig_reportfile.write(f'{log_entry}\n')
                    return True

        return False

    def _write_genes_with_stitched_contig(self, data):
        """
        Writes a file listing genes for which a stitched_contig was created (or skipped of flag --no_stitched_contig was
        provided to hybpiper assemble). These per-sample files are collated in the assemble.py script after
        all genes have completed.
        """

        with open(f'{self.prefix}/genes_with_stitched_contig.csv', 'w') as stitched_contig_reportfile:
            stitched_contig_reportfile.write(f'{data}\n')

    def _no_stitched_contig(self):
        """
        Identifies the single longest Exonerate hit; does not attempt to stitch multiple hits together into a
        stitched contig.

        :return Bio.SeqRecord.SeqRecord: SeqRecord for the single longest Exonerate hit.
        """

        if not self.hits_filtered_by_pct_similarity_dict:
            return None

        sample_name = os.path.split(self.prefix)[-1]
        gene_name = os.path.split(self.prefix)[-2]

        # Don't overwrite self.hits_subsumed_hits_removed_dict:
        exonerate_hits_subsumed_hits_removed_copy = copy.deepcopy(self.hits_subsumed_hits_removed_dict)

        sorted_by_hit_length = sorted(exonerate_hits_subsumed_hits_removed_copy.values(),
                                      key=lambda x: len(x['hit_sequence']), reverse=True)

        if len(sorted_by_hit_length) == 0:
            raise ValueError(f'The list sorted_by_hit_length for gene {gene_name} is empty!')

        sorted_by_hit_length[0]['hit_sequence'].description = f'Flag no_stitched_contig used. Single longest hit ' \
                                                              f'{sorted_by_hit_length[0]["hit_sequence"].id}'
        sorted_by_hit_length[0]['hit_sequence'].id = sample_name
        sorted_by_hit_length[0]['hit_sequence'].name = sample_name

        # Write report file:
        log_entry = f'{sample_name},{gene_name}, Stitched contig step skipped (user provided the ' \
                    f'--no_stitched_contig flag to hybpiper assemble). Gene sequence contains the longest Exonerate ' \
                    f'hit sequence only.'
        self._write_genes_with_stitched_contig(log_entry)

        return sorted_by_hit_length[0]['hit_sequence']

    def write_no_stitched_contig(self):
        """
        Writes DNA (suffix '.FNA') and amino-acid (suffix '.FAA') sequences for the single longest Exonerate hit (or
        single remaining hit sequence) to fasta file.

        :return NoneType: no explicit return
        """

        gene_name = os.path.split(self.prefix)[-2]
        hit_id = self.no_stitched_contig_seqrecord.id
        name = self.no_stitched_contig_seqrecord.name
        description = self.no_stitched_contig_seqrecord.description

        dna_seqrecord_to_write = self.no_stitched_contig_seqrecord
        amino_acid_seq_to_write = self.no_stitched_contig_seqrecord.seq.translate()
        amino_acid_seqrecord_to_write = SeqRecord(seq=amino_acid_seq_to_write, id=hit_id, name=name,
                                                  description=description)

        with open(f'{self.prefix}/sequences/FNA/{gene_name}.FNA', 'w') as fna_handle:
            SeqIO.write(dna_seqrecord_to_write, fna_handle, 'fasta')

        with open(f'{self.prefix}/sequences/FAA/{gene_name}.FAA', 'w') as faa_handle:
            SeqIO.write(amino_acid_seqrecord_to_write, faa_handle, 'fasta')

    def _group_hits_by_depth(self):  # try to split out contigs belongs to individual paralogs, rather than subsuming
        raise NotImplementedError

    def __repr__(self):
        """
        Returns a human-readable summary of the Exonerate object.
        """

        attrs = []
        for key in sorted(self.__dict__):
            attrs.append(f'{key}: {getattr(self, key)}')
        attrs_formatted = '\n'.join(attrs)
        return f'\nExonerate object:\n{"-" * 20}\n{attrs_formatted}\n{"-" * 20}\n'


def pairwise(iterable):
    """
    s -> (s0,s1), (s1,s2), (s2, s3), ...
    Used in the function fullContigs to iterate over overlapping pairs of hit_start_and_end_indices.
    """

    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def pairwise_longest(iterable):
    """
    s -> (s0,s1), (s1,s2), (s2, s3), ...
    Used in the function fullContigs to iterate over overlapping pairs of hit_start_and_end_indices.
    """

    a, b = tee(iterable)
    next(b, None)
    return itertools.zip_longest(a, b)


def grouped(iterable, n):
    """
    s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ...
    Used in the function fullContigs to iterate over non-overlapping pairs of reads from a sam file (i.e. reads 1+2,
    then reads 3+4 etc).
    """

    return zip(*[iter(iterable)] * n)


def set_stitched_contig_chimera_test(no_stitched_contig_bool, prefix):
    """
    Return True if a file of R1/R2 interleaved reads is found. Also return the path to the
    interleaved reads file.

    :param bool no_stitched_contig_bool: if True, no chimera test will be performed
    :param str prefix: path of gene/sample name
    :return: bool, str: path to interleaved fasta file for gene
    """

    logger = logging.getLogger(f'{os.path.split(prefix)[0]}')

    if not no_stitched_contig_bool:
        gene_folder = os.path.split(prefix)[0]
        interleaved_reads = f'{gene_folder}/{gene_folder}_interleaved.fasta'

        try:
            with open(interleaved_reads):
                pass
            return True, interleaved_reads
        except FileNotFoundError:
            logger.debug(f'Supercontigs will be generated, but there is no file of interleaved reads. Presuming '
                         f'single-end reads were provided: no chimera testing will be performed')
            return False, None
    else:
        return False, None


def parse_spades_and_best_reference(assemblyfile, proteinfile, prefix):
    """
    Return a SeqIO dictionary for the SPAdes contigs (assemblyfile) and the 'best' protein
    reference for the sample/gene (proteinfile).

    :param str assemblyfile: name of the FASTA file containing DNA sequence assembly
    :param str proteinfile: name of the FASTA file containing the target protein sequence
    :param str prefix: path of gene/sample name
    :return dict, dict spades_assembly_dict, best_protein_ref_dict:
    """

    logger = logging.getLogger(f'{os.path.split(prefix)[0]}')

    try:
        assemblyfile = open(assemblyfile)
    except IOError:
        logger.debug(f'The file {assemblyfile} could not be opened!')
        return
    try:
        proteinfile = open(proteinfile)
    except IOError:
        logger.debug(f'The file {proteinfile} could not be opened!')
        return

    spades_assembly_dict = SeqIO.to_dict(SeqIO.parse(assemblyfile, 'fasta'))
    best_protein_ref_dict = SeqIO.to_dict(SeqIO.parse(proteinfile, 'fasta'))

    return spades_assembly_dict, best_protein_ref_dict


def create_output_directories(prefix, assemblyfile):
    """

    :param str prefix: name of the sample
    :param str assemblyfile: name of the FASTA file containing DNA sequence assembly
    :return str: prefix: path of gene/sample name
    """

    if prefix:
        if os.path.exists(prefix):
            pass
        else:
            os.mkdir(prefix)
    else:
        prefix = os.path.basename(assemblyfile).split('.')[0]
        os.mkdir(prefix)

    if not os.path.exists(f'{prefix}/sequences/FNA'):
        os.makedirs(f'{prefix}/sequences/FNA')
    if not os.path.exists(f'{prefix}/sequences/FAA'):
        os.makedirs(f'{prefix}/sequences/FAA')

    return prefix


def standalone():
    """
    Used when this module is run as a stand-alone script. Parses command line arguments and runs function main().:
    """

    parser = argparse.ArgumentParser(
        description="exonerate_hits.py; Generate gene-by-gene protein and nucleotide files from Bait Capture Assembly")
    parser.add_argument("--debug", help="Print debugging information for development testing.",
                        action="store_true", dest="loglevel", default=False)
    parser.add_argument("proteinfile", help="FASTA file containing one 'target' sequence per protein.")
    parser.add_argument("assemblyfile", help="FASTA file containing DNA sequence assembly.")
    parser.add_argument("--prefix", help="Prefix for directory, files, and sequences generated from this assembly. If "
                                         "not specified, will be extracted from assembly file sample_name.",
                        default=None)
    parser.add_argument("--no_sequences", help="Do not generate protein and nucleotide sequence files.",
                        action="store_true", default=False)
    parser.add_argument("--first_search_filename",
                        help="Location of previously completed Exonerate results. Useful for testing.", default="no")
    parser.add_argument("-t", "--thresh",
                        help="Threshold for Percent Identity between contigs and proteins. default = 55%%", default=55,
                        type=int)
    parser.add_argument("--depth_multiplier",
                        help="Accept any full-length hit if it has a coverage depth X times the next best hit. Set to "
                             "zero to not use depth. Default = 10", default=10, type=int)
    parser.add_argument("--no_stitched_contig",
                        help="Do not create any stitched contigs. The longest single Exonerate hit will be used",
                        action="store_true", dest='no_stitched_contig', default=False)
    parser.add_argument("--bbmap_memory", help="memory (RAM ) to use for bbmap.sh", default=1, type=int)
    parser.add_argument("--bbmap_threads", help="threads to use for bbmap.sh", default=2, type=int)
    parser.add_argument("--bbmap_subfilter", default=7, type=int,
                        help="Ban bbmap.sh alignments with more than this many substitutions. Default is %(default)s")
    parser.add_argument("--chimeric_stitched_contig_edit_distance",
                        help="Minimum number of differences between one read of a read pair vs the stitched_contig "
                             "reference for a read pair to be flagged as discordant", default=7, type=int)
    parser.add_argument("--chimeric_stitched_contig_discordant_reads_cutoff",
                        help="minimum number of discordant reads pairs required to flag a stitched_contig as a "
                             "potential chimera of contigs from multiple paralogs", default=100, type=int)
    parser.add_argument("--paralog_warning_min_length_percentage", default=0.75, type=float,
                        help="Minimum length percentage of a contig vs reference protein length for a paralog warning "
                             "to be generated. Default is %(default)s")
    parser.add_argument("--pad_stitched_contig_query_gaps_with_n",
                        help="When constructing stitched contigs, pad any gaps between hits (with respect to the "
                             "query protein) with a number of Ns corresponding to the query gap multiplied by 3",
                        action="store_true", dest='stitched_contig_pad_n', default=False)
    parser.add_argument("--run_intronerate",
                        help="Run intronerate to recover fasta files for supercontig with introns (if present), "
                             "and introns-only", action="store_true", dest='intronerate', default=False)
    parser.add_argument("--no_padding_supercontigs",
                        help='If Intronerate is run, and a supercontig is created by concatenating multiple SPAdes '
                             'contigs, do not add 10 "N" characters between contig joins. By default, Ns will be '
                             'added.', action='store_true', dest='no_padding_supercontigs',
                        default=False)
    parser.add_argument('--keep_intermediate_files',
                        help='Keep all intermediate files and logs, which can be useful for '
                             'debugging. Default action is to delete them, which greatly reduces the total file '
                             'number).',
                        action='store_true', dest='keep_intermediate_files', default=False)

    args = parser.parse_args()

    main(args)


########################################################################################################################
# Define main()
########################################################################################################################

def main(args):

    # Setup logger for main():
    logger = logging.getLogger()
    ch = logging.StreamHandler()
    logger.addHandler(ch)
    if args.loglevel:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)

    # Create directories for output files based on the prefix name, or assemblyfile name:
    prefix = create_output_directories(args.prefix, args.assemblyfile)

    # Set whether the chimeric stitched contig test will be performed, and whether a file of interleaved reads is found:
    perform_stitched_contig_chimera_test, path_to_interleaved_fasta = set_stitched_contig_chimera_test(
        args.no_stitched_contig,
        prefix)

    # Read the SPAdes contigs and the 'best' protein reference seq into SeqIO dictionaries:
    spades_assembly_dict, best_protein_ref_dict = parse_spades_and_best_reference(args.assemblyfile,
                                                                                  args.proteinfile,
                                                                                  prefix)

    # Perform Exonerate search with 'best' protein ref as query and SPAdes contigs as subjects
    exonerate_text_output = initial_exonerate(args.proteinfile,
                                              args.assemblyfile,
                                              prefix)

    exonerate_result = parse_exonerate_and_get_stitched_contig(exonerate_text_output,
                                                               query_file=args.proteinfile,
                                                               paralog_warning_min_length_percentage=
                                                               args.paralog_warning_min_length_percentage,
                                                               thresh=args.thresh,
                                                               logger=logger,
                                                               prefix=prefix,
                                                               discordant_cutoff=
                                                               args.chimeric_stitched_contig_discordant_reads_cutoff,
                                                               edit_distance=args.chimeric_stitched_contig_edit_distance,
                                                               bbmap_subfilter=args.bbmap_subfilter,
                                                               bbmap_memory=args.bbmap_memory,
                                                               bbmap_threads=args.bbmap_threads,
                                                               interleaved_fasta_file=path_to_interleaved_fasta,
                                                               no_stitched_contig=args.no_stitched_contig,
                                                               spades_assembly_dict=spades_assembly_dict,
                                                               depth_multiplier=args.depth_multiplier,
                                                               keep_intermediate_files=args.keep_intermediate_files)
    if not exonerate_result.stitched_contig_seqrecord:
        return

    logger.debug(f'There were {len(exonerate_result.hits_filtered_by_pct_similarity_dict)} Exonerate '
                 f'hits for {args.proteinfile} after filtering by similarity threshold {args.thresh}.')

    if intronerate and exonerate_result and exonerate_result.hits_filtered_by_pct_similarity_dict:
        logger.debug(f'exonerate_result.hits_subsumed_hits_removed_overlaps_trimmed_dict is:'
                     f' {exonerate_result.hits_subsumed_hits_removed_overlaps_trimmed_dict}')
        logger.debug(f'Running intronerate')
        intronerate(exonerate_result,
                    spades_assembly_dict,
                    logger=logger,
                    no_padding_supercontigs=args.no_padding_supercontigs,
                    keep_intermediate_files=args.keep_intermediate_files)


########################################################################################################################
# Run the script
########################################################################################################################
if __name__ == "__main__":
    standalone()

################################################## END OF SCRIPT #######################################################
