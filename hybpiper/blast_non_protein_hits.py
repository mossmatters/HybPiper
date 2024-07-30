#!/usr/bin/env python

"""
This module takes a nucleotide SPAdes assembly for a non-protein-coding locus and a nucleotide sequence for the
'best' target file sequence for the locus. It performs a BLASTn search using the targetfile nucleotide sequence as the
query and the SPAdes contigs as targets and extracts hit regions to create (where possible) a contiguous coding
sequence for the gene/sample.

The module also performs the following tasks:

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
import textwrap
import numpy as np


def initial_blast(locusfilename,
                  assemblyfilename,
                  prefix,
                  blast_contigs_task,
                  keep_intermediate_files=False):
    """
    Conduct BLASTn search.

    :param str locusfilename: path to the chosen target-file locus query fasta file
    :param str assemblyfilename: path to the SPAdes assembly contigs file
    :param str prefix: path of gene/sample name
    :param str blast_contigs_task: task when running blastn searches
    :param bool keep_intermediate_files: if True, keep intermediate files
    :return NoneType or str, None or outputfilename:
        - None: returned if BLASTn searches fail to produce output file
        - outputfilename: the BLASTn xml file output
    """

    logger = logging.getLogger(f'{os.path.split(prefix)[0]}')

    blast_database_name = os.path.splitext(os.path.basename(assemblyfilename))[0]
    blast_database_full_path = f'{prefix}/{blast_database_name}'

    # Create a BLAST database for the SPAdes assemblies:
    makeblastdb_cmd = f'makeblastdb -dbtype nucl -in {assemblyfilename} -out {blast_database_full_path}'

    try:
        result = subprocess.run(makeblastdb_cmd,
                                shell=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                universal_newlines=True,
                                check=True)

        logger.debug(f'makeblastdb check_returncode() is: {result.check_returncode()}')
        logger.debug(f'makeblastdb stdout is: {result.stdout}')
        logger.debug(f'makeblastdb stderr is: {result.stderr}')

    except subprocess.CalledProcessError as exc:
        logger.error(f'makeblastdb FAILED. Output is: {exc}')
        logger.error(f'makeblastdb stdout is: {exc.stdout}')
        logger.error(f'makeblastdb stderr is: {exc.stderr}')
        return None

    # Run the BLASTn search:
    outputfilename = f'{prefix}/blastn_results.xml'

    blastn_command = (f'blastn '
                      f'-task {blast_contigs_task} '
                      f'-db {blast_database_full_path} '
                      f'-query {locusfilename} '
                      f'-outfmt 5 '
                      f'-out {outputfilename}')

    logger.debug(f'BLASTn command is: {blastn_command}')

    try:
        subprocess.run(blastn_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       universal_newlines=True, check=True)
    except subprocess.CalledProcessError as exc:
        logger.debug(f'BLASTn FAILED for {prefix}. Output is: {exc}')
        logger.debug(f'BLASTn stdout is: {exc.stdout}')
        logger.debug(f'BLASTn stderr is: {exc.stderr}')

    # Delete blast db files:
    blastdb_files = [
        # 'f'{prefix}/{blast_database_name}.'
    ]



    if utils.file_exists_and_not_empty(outputfilename):
        logger.debug('BLASTn ran successfully')
        return outputfilename
    else:
        return None


def parse_blast_and_get_stitched_contig(blast_xml_output,
                                        query_file,
                                        paralog_warning_min_length_percentage,
                                        thresh,
                                        logger,
                                        prefix,
                                        chimera_check,
                                        discordant_cutoff,
                                        edit_distance,
                                        bbmap_subfilter,
                                        bbmap_memory,
                                        bbmap_threads,
                                        interleaved_fasta_file,
                                        no_stitched_contig,
                                        stitched_contig_pad_n,
                                        spades_assembly_dict,
                                        depth_multiplier,
                                        keep_intermediate_files,
                                        blast_hit_sliding_window_size,
                                        blast_hit_sliding_window_thresh,
                                        verbose_logging):
    """
    => Parses the C4 alignment text output of Exonerate using BioPython SearchIO.
    => Generates paralog warning and fasta files.
    => Generates stitched contig (or single hit) fasta files for nucleotide and amino-acid sequences.
    => Performs a stitched contig chimera test if file of R1/R2 interleaved reads is present and no_stitched_contig is
    False.

    :param str blast_xml_output: path to the results xml file output by BLASTn
    :param str query_file: path to the protein query fasta file.
    :param float paralog_warning_min_length_percentage: percentage coverage of query required for paralog warning.
    :param int thresh: minimum percentage similarity threshold used to filter BLASTn hits.
    :param logger logger: a logger object
    :param str prefix: path to gene/sample folder e.g. gene001/sampleID
    :param bool chimera_check: run chimera check. Default is False
    :param int discordant_cutoff: number of discordant read pairs for a stitched contig to be flagged as chimeric
    :param int edit_distance: edit distance threshold for identifying discordant read pairs
    :param int bbmap_subfilter: ban bbmap.sh alignments with more than this many substitutions
    :param int bbmap_memory: MB of RAM to use for bbmap.sh
    :param int bbmap_threads: number of threads to use for bbmap.sh
    :param None, str interleaved_fasta_file: path to the file of interleaved R1 and R2 fasta seqs, if present
    :param bool no_stitched_contig: if True, return the longest BLASTn hit only
    :param bool stitched_contig_pad_n: if True, pad gaps in stitched contig with Ns corresponding to query gap * 3
    :param dict spades_assembly_dict: a dictionary of raw SPAdes contigs
    :param int depth_multiplier: assign long paralog as main if coverage depth <depth_multiplier> time other paralogs
    :param bool keep_intermediate_files: if True, keep intermediate files from stitched contig
    :param int blast_hit_sliding_window_size: size of the sliding window (in nucletides) when trimming termini
    of BLASTn hits
    :param int blast_hit_sliding_window_thresh: percentage similarity threshold for the sliding window (in
    nucleotides) when trimming termini of BLASTn hits
    :param bool verbose_logging: if True, log additional information to file

    :return __main__.BLAST_CONTIGS: instance of the class BLAST_CONTIGS for a given gene
    """

    try:
        blast_hits_from_alignment = list(SearchIO.parse(blast_xml_output, 'blast-xml'))
    except AttributeError as e:
        logger.debug(f'Error parsing BLASTn results file for prefix {prefix}\nError is: {e} ')
        return None

    logger.debug(f'no_stitched_contig is: {no_stitched_contig}')

    blast_result = BlastContigs(searchio_object=blast_hits_from_alignment,
                                query_file=query_file,
                                paralog_warning_min_length_percentage=paralog_warning_min_length_percentage,
                                thresh=thresh,
                                logger=logger,
                                prefix=prefix,
                                chimera_check=chimera_check,
                                discordant_cutoff=discordant_cutoff,
                                edit_distance=edit_distance,
                                bbmap_subfilter=bbmap_subfilter,
                                bbmap_memory=bbmap_memory,
                                bbmap_threads=bbmap_threads,
                                interleaved_fasta_file=interleaved_fasta_file,
                                no_stitched_contig=no_stitched_contig,
                                stitched_contig_pad_n=stitched_contig_pad_n,
                                spades_assembly_dict=spades_assembly_dict,
                                depth_multiplier=depth_multiplier,
                                keep_intermediate_files=keep_intermediate_files,
                                blast_hit_sliding_window_size=blast_hit_sliding_window_size,
                                blast_hit_sliding_window_thresh=blast_hit_sliding_window_thresh,
                                verbose_logging=verbose_logging)

    if verbose_logging:
        logger.debug(blast_result)

    blast_result.write_exonerate_stats_file()

    if not blast_result.hits_subsumed_hits_removed_dict:  # i.e. no hits left after filtering
        return None

    # if blast_result.long_paralogs_dict:  # i.e. there are long paralogs recovered
    blast_result.write_long_paralogs_and_warnings_to_file()

    if no_stitched_contig:
        blast_result.write_no_stitched_contig()
    else:
        blast_result.write_stitched_contig_to_file()

    if keep_intermediate_files:
        blast_result.write_trimmed_stitched_contig_hits_to_file()

    return blast_result


class BlastContigs(object):
    """
    Class to parse BLASTn XML results (SearchIO object) for a given gene. Returns a Blast_Contigs object.
    """

    def __init__(self,
                 searchio_object,
                 query_file=None,
                 paralog_warning_min_length_percentage=0.75,
                 thresh=55,
                 logger=None,
                 prefix=None,
                 chimera_check=False,
                 discordant_cutoff=5,
                 edit_distance=5,
                 bbmap_subfilter=7,
                 bbmap_memory=1,
                 bbmap_threads=1,
                 interleaved_fasta_file=None,
                 no_stitched_contig=False,
                 stitched_contig_pad_n=True,
                 spades_assembly_dict=None,
                 depth_multiplier=10,
                 keep_intermediate_files=False,
                 blast_hit_sliding_window_size=9,
                 blast_hit_sliding_window_thresh=65,
                 verbose_logging=False):
        """
        Initialises class attributes.

        :param list searchio_object: list returned by parsing BLASTn output with SearchIO.parse
        :param str query_file: path to the nucleotide query fasta file
        :param float paralog_warning_min_length_percentage: percentage coverage of query required for paralog warning
        :param int thresh: minimum percentage similarity threshold used to filter BLASTn hits
        :param logging.RootLogger logger: a logger object
        :param str prefix: path to gene/sample folder e.g. gene001/sampleID
        :param bool chimera_check: run chimera check. Default is False
        :param int discordant_cutoff: number of discordant read pairs for a stitched contig to be flagged as chimeric
        :param int edit_distance: edit distance threshold for identifying discordant read pairs
        :param int bbmap_subfilter: ban bbmap.sh alignments with more than this many substitutions
        :param int bbmap_memory: MB of RAM to use for bbmap.sh
        :param int bbmap_threads: number of threads to use for bbmap.sh
        :param str interleaved_fasta_file: path to the file of interleaved R1 and R2 fasta seqs, if present
        :param bool no_stitched_contig: if True, return the longest Exonerate hit only
        :param bool stitched_contig_pad_n: if True, pad gaps in stitched contig with Ns corresponding to query gap * 3
        :param dict spades_assembly_dict: a dictionary of raw SPAdes contigs
        :param int depth_multiplier: assign long paralog as main if coverage depth <depth_multiplier> other paralogs
        :param bool keep_intermediate_files: if True, keep intermediate files from stitched contig
        :param int blast_hit_sliding_window_size: size of the sliding window (in nucleotides) when trimming
        termini of BLASTn hits
        :param int blast_hit_sliding_window_thresh: percentage similarity threshold for the sliding window (
        in nucleotides) when trimming termini of BLASTn hits
        :param bool verbose_logging: if True, log additional information to file
        """

        if len(searchio_object) != 1:  # This should always be 1 for a single BLASTn query
            raise ValueError(f'searchio_object list is greater than 1!')

        self.blast_searchio_alignment = searchio_object
        self.query_id = searchio_object[0].id
        self.query_length = len(SeqIO.read(query_file, 'fasta'))
        self.similarity_threshold = thresh
        self.sliding_window_size = blast_hit_sliding_window_size
        self.sliding_window_thresh = blast_hit_sliding_window_thresh
        self.paralog_warning_by_contig_length_pct = paralog_warning_min_length_percentage
        self.logger = logger
        self.prefix = prefix
        self.no_stitched_contig = no_stitched_contig
        self.stitched_contig_pad_n = stitched_contig_pad_n
        self.interleaved_fasta_file = interleaved_fasta_file
        self.chimera_discordant_cutoff = discordant_cutoff
        self.chimera_edit_distance = edit_distance
        self.chimera_bbmap_subfilter = bbmap_subfilter
        self.chimera_bbmap_memory = bbmap_memory
        self.chimera_bbmap_threads = bbmap_threads
        self.spades_assembly_dict = spades_assembly_dict
        self.depth_multiplier = depth_multiplier
        self.keep_intermediate_files = keep_intermediate_files
        self.verbose_logging = verbose_logging
        self.hits_filtered_by_pct_similarity_dict = self._parse_searchio_object()
        self.hits_subsumed_hits_removed_dict = self._remove_subsumed_hits()
        self.hits_subsumed_hits_removed_overlaps_trimmed_dict = self._trim_overlapping_hits()
        self.long_paralogs_dict = self._recover_long_paralogs()
        self.paralog_warning_by_contig_depth = self._paralog_warning_by_contig_depth()
        self.stitched_contig_seqrecord = self._create_stitched_contig()
        self.stitched_contig_hit_ranges = self._get_stitched_contig_hit_ranges()
        self.no_stitched_contig_seqrecord = self._no_stitched_contig() if self.no_stitched_contig else None

        # If chimera_check is True, only perform chimera test if stitched contigs are being created AND
        # interleaved_fasta_file is not None AND a multi-hit stitched contig has been created:
        if chimera_check:
            if (self.hits_filtered_by_pct_similarity_dict
                    and not self.no_stitched_contig
                    and interleaved_fasta_file
                    and not self.stitched_contig_seqrecord.description == 'single_hit'):
                self.chimera_warning_bool = self._stitched_contig_chimera_warning()
            else:
                self.chimera_warning_bool = None
        else:
            self.chimera_warning_bool = None

    def _parse_searchio_object(self):
        """
        Parses the object returned by BioPython SearchIO.parse.
        => Calculates query-vs-hit similarity scores for each hit, and filters hit based on a given threshold
        => Sorts similarity-filtered hsps by start position then end position in the protein query
        => Applies a sliding window similarity filter to 5' and 3' hit termini, and trims if below similarity threshold
        => Populates a dict of dicts for each hit, with hitname: {key:value hit data}

        :return collections.defaultdict filtered_by_similarity_hsps_dict: dict of dicts for each hit
        """

        single_exonerate_qresult = self.blast_searchio_alignment[0]
        filtered_hsps = []

        ############################################################################################################
        # Calculate hsp similarity and filter hsps via a given threshold similarity percentage:
        ############################################################################################################
        for hsp in single_exonerate_qresult.hsps:

            similarity_count_total = 0
            similarity_count = 0

            for alignment in hsp.aln_annotation_all:  # Note that for SearchIO blast xml, only contains similarity
                for alignment_column in alignment['similarity']:
                    if alignment_column == '|':
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

        ############################################################################################################
        # Sort hsps list by query start location THEN query end location:
        ############################################################################################################
        for filtered_hsp in sorted(filtered_hsps, key=lambda x: (x[0].query_start, x[0].query_end)):

            # Trim ends of filtered-for-similarity hsps if they fall beneath self.similarity_threshold. First, check if
            # there is more than one fragment in the hsp. If so, recover and concatenate the similarity annotations
            # and hit (nucleotide codon) annotations:
            hsp = filtered_hsp[0]
            assert len(hsp.aln_annotation_all) == 1

            hit_similarity_symbols = hsp.aln_annotation_all[0]['similarity']
            hit_bases = hsp.hit.seq

            ############################################################################################################
            # Calculate similarity values within a sliding window:
            ############################################################################################################
            window_size = self.sliding_window_size
            window_thresh = self.sliding_window_thresh
            window_similarity_percentages = []
            all_window_similarity_symbols = []
            all_window_bases = []

            for i in range(0, len(hit_similarity_symbols) - (window_size - 1)):
                window_similarity_symbols = hit_similarity_symbols[i:i + window_size]
                window_bases = hit_bases[i:i + window_size]
                all_window_similarity_symbols.append(window_similarity_symbols)
                all_window_bases.append(window_bases)
                window_identical_bases_count = 0
                window_all_bases_count = 0

                for symbol in window_similarity_symbols:
                    window_all_bases_count += 1
                    if symbol == '|':
                        window_identical_bases_count += 1

                try:
                    window_similarity = \
                        f'{window_identical_bases_count / window_all_bases_count * 100:.2f}'
                except ZeroDivisionError:
                    window_similarity = float(0.0)

                window_similarity_percentages.append(float(window_similarity))

            ############################################################################################################
            # Calculate indices of where a sliding-window similarity value crosses the similarity threshold:
            ############################################################################################################
            data = np.array(window_similarity_percentages)
            data_adjusted = np.where(data <= window_thresh, 0, 1)  # 0 or 1 if above or below/equal

            # Below: prepend=1 means that the start value (before first sliding window value) corresponds to good
            # quality sequence. So, if the first sliding window value is good, a value of 0 (no change) will be
            # returned, and no upward crossing will be detected:
            five_prime_threshold_crossings = np.diff(data_adjusted >= 1, prepend=1)
            three_prime_threshold_crossings = np.diff(data_adjusted[::-1] >= 1, prepend=1)  # reverse list)

            # Below: every second item starting at second value, from first index (0) in each element:
            five_prime_upward_crossings = np.argwhere(five_prime_threshold_crossings)[1::2, 0]
            # Below: every second item starting at first value, from first index (0) in each element:
            five_prime_downward_crossings = np.argwhere(five_prime_threshold_crossings)[::2, 0]
            three_prime_upward_crossings = np.argwhere(three_prime_threshold_crossings)[1::2, 0]
            three_prime_downward_crossings = np.argwhere(three_prime_threshold_crossings)[::2, 0]

            if self.verbose_logging:
                self.logger.debug(f'five_prime_upward_crossings: {five_prime_upward_crossings}')
                self.logger.debug(f'five_prime_downward_crossings: {five_prime_downward_crossings}')
                self.logger.debug(f'three_prime_upward_crossings: {three_prime_upward_crossings}')
                self.logger.debug(f'three_prime_downward_crossings: {three_prime_downward_crossings}')

            ############################################################################################################
            # Recover the first upwards and downwards crossing:
            ############################################################################################################

            # Five prime:
            try:
                five_prime_first_upward_crossing_nucleotides = five_prime_upward_crossings[0]
            except IndexError:
                five_prime_first_upward_crossing_nucleotides = 'no_crossing'
            try:
                five_prime_first_downward_crossing_nucleotides = five_prime_downward_crossings[0]
            except IndexError:
                five_prime_first_downward_crossing_nucleotides = 'no_crossing'

            # Three prime:
            try:
                three_prime_first_upward_crossing_nucleotides = three_prime_upward_crossings[0]
            except IndexError:
                three_prime_first_upward_crossing_nucleotides = 'no_crossing'
            try:
                three_prime_first_downward_crossing_nucleotides = three_prime_downward_crossings[0]
            except IndexError:
                three_prime_first_downward_crossing_nucleotides = 'no_crossing'

            ############################################################################################################
            # Get slice indices, if any:
            ############################################################################################################
            five_prime_slice = 0  # i.e. default is the first position in the sequence
            three_prime_slice = 0  # i.e. default is the first (converted from last) position in the sequence

            # Check if 5' trimming is required:
            if five_prime_first_downward_crossing_nucleotides == 'no_crossing':  # no sliding window beneath the
                # self.similarity_threshold value
                if self.verbose_logging:
                    self.logger.debug(f'five_prime_first_downward_crossing_nucleotides is'
                                      f' {five_prime_first_downward_crossing_nucleotides}: no trimming required for '
                                      f'prefix {self.prefix} hsp {hsp.hit_id}')

            elif five_prime_first_downward_crossing_nucleotides != 0:  # i.e. hit starts above threshold
                if self.verbose_logging:
                    self.logger.debug(f'five_prime_first_downward_crossing_nucleotides is'
                                      f' {five_prime_first_downward_crossing_nucleotides}; hit starts above '
                                      f'threshold, so no 5-prime trimming required for prefix {self.prefix} hsp '
                                      f'{hsp.hit_id}')
            else:
                five_prime_slice = five_prime_first_upward_crossing_nucleotides
                if self.verbose_logging:
                    self.logger.debug(f'5-prime trimming for prefix {self.prefix} hsp {hsp.hit_id} is '
                                      f'{five_prime_slice}')

            # Check if 3' trimming is required:
            if three_prime_first_downward_crossing_nucleotides == 'no_crossing':
                if self.verbose_logging:
                    self.logger.debug(f'three_prime_first_downward_crossing_nucleotides is '
                                      f'{three_prime_first_downward_crossing_nucleotides}: no trimming required for '
                                      f'prefix {self.prefix} hsp {hsp.hit_id}')

            elif three_prime_first_downward_crossing_nucleotides != 0:  # i.e. hit starts above threshold
                if self.verbose_logging:
                    self.logger.debug(f'three_prime_first_downward_crossing_nucleotides is'
                                      f' {three_prime_first_downward_crossing_nucleotides}; hit starts above '
                                      f'threshold, so no 3-prime trimming required for prefix {self.prefix} hsp '
                                      f'{hsp.hit_id}')
            else:
                three_prime_slice = three_prime_first_upward_crossing_nucleotides
                if self.verbose_logging:
                    self.logger.debug(f'3-prime trimming for prefix {self.prefix} hsp {hsp.hit_id} is '
                                      f'{three_prime_slice}')

            ############################################################################################################
            # Adjust ends of slices to start with the first nucleotide with similarity '|':
            ############################################################################################################

            if five_prime_slice:  # i.e., five_prime_slice is not zero
                window_similarity_five_prime_slice = all_window_similarity_symbols[five_prime_slice]
                for symbol in window_similarity_five_prime_slice:
                    if symbol != '|':
                        five_prime_slice = five_prime_slice + 1
                    else:
                        break

            if three_prime_slice:  # i.e., it's not zero
                window_similarity_three_prime_slice = \
                    all_window_similarity_symbols[::-1][three_prime_slice][::-1]
                for symbol in window_similarity_three_prime_slice:
                    if symbol != '|':
                        three_prime_slice = three_prime_slice + 1
                    else:
                        break

            ############################################################################################################
            # Re-calculate the hit similarity based on the sliced sequence:
            ############################################################################################################
            hit_similarity_slice = \
                hit_similarity_symbols[five_prime_slice: len(hit_similarity_symbols) - three_prime_slice]

            similarity_count_total = 0
            similarity_count = 0

            for symbol in hit_similarity_slice:
                if symbol == '|':
                    similarity_count += 1
                similarity_count_total += 1
            hit_similarity = f'{similarity_count / similarity_count_total * 100:.2f}'

            ############################################################################################################
            # Extract values from hsp object and adjust as necessary for trim slices:
            ############################################################################################################
            spades_contig_depth = float(filtered_hsp[0].hit_id.split('_')[-1])  # dependant on SPAdes header output!
            query_range_original = filtered_hsp[0].query_range
            query_range = (round(query_range_original[0] + five_prime_slice),
                           round((query_range_original[1] - three_prime_slice)))
            query_range_all = filtered_hsp[0].query_range_all
            hit_inter_ranges = filtered_hsp[0].hit_inter_ranges
            hit_inter_spans = filtered_hsp[0].hit_inter_spans
            hsp_hit_strand_all = filtered_hsp[0].hit_strand_all
            assert len(set(hsp_hit_strand_all)) == 1  # Check that all HSP fragments are on the same strand
            hsp_hit_strand = next(iter(set(hsp_hit_strand_all)))
            hit_range_original = filtered_hsp[0].hit_range
            hit_range = (round(hit_range_original[0] + five_prime_slice),
                         round((hit_range_original[1] - three_prime_slice)))
            hit_range_all_original = filtered_hsp[0].hit_range_all

            # Adjust hit_range_all_original for 5' and 3' slices, if any:
            hit_range_all_flattened = [item for sublist in hit_range_all_original for item in sublist]
            if hsp_hit_strand == -1:
                if len(hit_range_all_original) == 1:  # i.e. a single coordinate tuple
                    hit_range_all_flattened[0] = round(hit_range_original[0] + five_prime_slice)
                    hit_range_all_flattened[-1] = round(hit_range_original[1] - three_prime_slice)
                else:  # more than one coordinate tuple
                    hit_range_all_flattened[-2] = round(hit_range_original[0] + five_prime_slice)
                    hit_range_all_flattened[1] = round(hit_range_original[1] - three_prime_slice)
            else:
                hit_range_all_flattened[0] = round(hit_range_original[0] + five_prime_slice)
                hit_range_all_flattened[-1] = round(hit_range_original[1] - three_prime_slice)

            hit_range_all = []  # re-make the hit_range_all list with adjusted values
            for pair in grouped(hit_range_all_flattened, 2):
                hit_range_all.append(pair)

            ############################################################################################################
            # Set a unique hit name for cases where there's >1 hit for a single SPAdes contig:
            ############################################################################################################
            unique_hit_name = f'{filtered_hsp[0].hit_id},{self.query_id},{query_range[0]},{query_range[1]}' \
                              f',{hit_similarity},({hsp_hit_strand}),{hit_range[0]},{hit_range[1]}'

            seq = Seq(hit_bases)  # Do _not_ strip gaps

            # Trim/slice the seq if required:
            if not three_prime_slice:
                seq_sliced = seq[five_prime_slice:]
            else:
                seq_sliced = seq[five_prime_slice: len(seq) - three_prime_slice]

            hit_seq = SeqRecord(id=unique_hit_name, name=unique_hit_name, description=unique_hit_name, seq=seq_sliced)

            # Populate nested dictionary for hsp:
            filtered_by_similarity_hsps_dict[unique_hit_name]['query_range_original'] = query_range_original
            filtered_by_similarity_hsps_dict[unique_hit_name]['query_range'] = query_range
            filtered_by_similarity_hsps_dict[unique_hit_name]['query_range_all'] = query_range_all
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_range_original'] = hit_range_original
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_range'] = hit_range
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_range_all_original'] = hit_range_all_original
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_range_all'] = hit_range_all
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_inter_ranges'] = hit_inter_ranges
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_inter_spans'] = hit_inter_spans
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_strand'] = hsp_hit_strand
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_sequence'] = hit_seq
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_spades_contig_depth'] = spades_contig_depth
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_similarity_original'] = filtered_hsp[1]
            filtered_by_similarity_hsps_dict[unique_hit_name]['hit_similarity'] = hit_similarity

        # The sliding window trim filter can change the query start/end order of hits, so re-sort the dictionary:
        filtered_by_similarity_hsps_dict_sorted = \
            {sorted_key: filtered_by_similarity_hsps_dict[sorted_key] for sorted_key in
             sorted(filtered_by_similarity_hsps_dict,
                    key=lambda x: (filtered_by_similarity_hsps_dict[x]['query_range'][0],
                                   filtered_by_similarity_hsps_dict[x]['query_range'][1]))}

        return filtered_by_similarity_hsps_dict_sorted

    def _remove_subsumed_hits(self):
        """
        Takes a dictionary of hits that have been filtered via similarity to the query, and removes any hit that has
        a query range that completely subsumes (i.e. encompasses) another hit. If two hits have an identical query
        range (and are not themselves subsumed by another longer hit) the hit with the highest similarity is retained.

        :return collections.defaultdict: exonerate_hits_filtered_no_subsumed
        """

        if not self.hits_filtered_by_pct_similarity_dict:
            return None

        blast_hits_filtered_no_subsumed = copy.deepcopy(self.hits_filtered_by_pct_similarity_dict)
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
                        if self.verbose_logging:
                            self.logger.debug(f'to_remove {to_remove} is already in seqs_removed: {seqs_removed}')
                        continue
                    seqs_removed.append(to_remove)

            # If hit not already removed, remove it from the dict:
            if to_remove:
                try:
                    del blast_hits_filtered_no_subsumed[to_remove]
                except KeyError:
                    if self.verbose_logging:
                        self.logger.debug(f'hit {to_remove} already removed from dict')
                    pass

        # Select a single sequence from each range in hits_with_identical_range_and_similarity_dict:
        if len(hits_with_identical_range_and_similarity_dict) != 0:
            self.logger.debug(f'Gene has hits with identical query ranges and similarities; selecting one hit for '
                              f'each range')
            if self.verbose_logging:
                self.logger.debug(f'Dictionary hits_with_identical_range_and_similarity_dict is:'
                                  f' {hits_with_identical_range_and_similarity_dict}')

            for query_range, hits in hits_with_identical_range_and_similarity_dict.items():
                to_remove = list(hits)[1:]  # arbitrarily retain first hit if range and similarity are the same
                for hit_to_remove in to_remove:
                    try:
                        del blast_hits_filtered_no_subsumed[hit_to_remove]
                    except KeyError:
                        if self.verbose_logging:
                            self.logger.debug(f'hit {hit_to_remove} already removed from dict')
                        pass

        return blast_hits_filtered_no_subsumed

    def _trim_overlapping_hits(self):

        """
        For constructing the coding-seq-only stitched contig via _create_stitched_contig()

        => Takes a dictionary of hits that has been filtered via hit similarity to the query, and has had subsumed hits
        removed. If any of the remaining hits have overlaps in query ranges, the 3' end of the left hit is trimmed to
        remove overlap sequence.
        => If there are any gaps between hit pairs with respect to the query sequence, pads the 3' end of the left hit
        with a corresponding number of ends Ns.
        => For any trimmed hit, annotate the description of the corresponding SeqRecord with '3prime overlap trimmed by:
        <int>'.

        :return collections.defaultdict: exonerate_hits_subsumed_and_trimmed_dict
        """

        if not self.hits_subsumed_hits_removed_dict:
            return None

        if len(self.hits_subsumed_hits_removed_dict) == 1:  # i.e. single hit remaining from previous filtering
            for key, value in self.hits_subsumed_hits_removed_dict.items():
                value['hit_sequence'].description = f'Single hit after filtering: N/A'
                seq_no_gaps = value['hit_sequence'].seq.replace('-', '')
                value['hit_sequence'].seq = seq_no_gaps
            return self.hits_subsumed_hits_removed_dict

        # Don't overwrite self.hits_subsumed_hits_removed_dict:
        blast_hits_subsumed_hits_removed_copy = copy.deepcopy(self.hits_subsumed_hits_removed_dict)

        blast_hits_subsumed_and_trimmed_dict = {}
        for pair in list(pairwise_longest(blast_hits_subsumed_hits_removed_copy.values())):
            if pair[1] is not None:  # as pairwise_longest() will pad a hit without a pair with 'None'
                left_seq_hit_name, right_seq_hit_name = pair[0]['hit_sequence'].id, pair[1]['hit_sequence'].id
                left_seq_query_range, right_seq_query_range = pair[0]['query_range'], pair[1]['query_range']
                left_seq_seq, right_seq_seq = pair[0]['hit_sequence'], pair[1]['hit_sequence']

                # If overlapping hits, always trim the 3' end of the left hit:
                if left_seq_query_range[1] > right_seq_query_range[0]:
                    num_nucleotides_overlap = left_seq_query_range[1] - right_seq_query_range[0]

                    seq_to_keep = left_seq_seq[: -num_nucleotides_overlap]
                    seq_to_trim = left_seq_seq[-num_nucleotides_overlap:]

                    assert len(seq_to_keep) + len(seq_to_trim) == len(left_seq_seq)

                    assert len(seq_to_keep) >= 3, (f'{self.prefix} left_seq_hit_name: {left_seq_hit_name}, '
                                                   f'right_seq_hit_name: {right_seq_hit_name}')

                    seq_to_keep_no_gaps = seq_to_keep.seq.replace('-', '')
                    seq_to_trim_no_gaps = seq_to_trim.seq.replace('-', '')

                    pair[0]['hit_sequence'].seq = Seq(seq_to_keep_no_gaps)
                    pair[0]['hit_sequence'].description = f'3prime overlap trimmed by: {len(seq_to_trim_no_gaps)}'
                    blast_hits_subsumed_and_trimmed_dict[left_seq_hit_name] = pair[0]

                else:  # if no overlap, add left hit unmodified (3' padded with Ns if gap before next hit):
                    bases_in_gap_between_hits = right_seq_query_range[0] - left_seq_query_range[1]
                    seq_no_gaps = pair[0]['hit_sequence'].seq.replace('-', '')

                    if self.stitched_contig_pad_n:
                        pair[0]['hit_sequence'].seq = Seq(f"{seq_no_gaps}{'N' * bases_in_gap_between_hits}")
                    else:
                        pair[0]['hit_sequence'].seq = seq_no_gaps

                    pair[0]['hit_sequence'].description = f'No overlap: N/A'
                    blast_hits_subsumed_and_trimmed_dict[left_seq_hit_name] = pair[0]

            else:  # process final unpaired left hit
                left_seq_hit_name = pair[0]['hit_sequence'].id
                pair[0]['hit_sequence'].description = f'No overlap: N/A'
                seq_no_gaps = pair[0]['hit_sequence'].seq.replace('-', '')
                pair[0]['hit_sequence'].seq = seq_no_gaps
                blast_hits_subsumed_and_trimmed_dict[left_seq_hit_name] = pair[0]

        return blast_hits_subsumed_and_trimmed_dict

    def write_trimmed_stitched_contig_hits_to_file(self):
        """
        Writes DNA (suffix '.FNA') sequences for the filtered, sorted and trimmed BLASTn hits to fasta file. Used
        for debugging.

        :return NoneType: no explicit return
        """

        # Write DNA seqs
        with open(f'{self.prefix}/blast_hits_trimmed.FNA', 'w') as fasta_handle_nucl:
            for key, value in self.hits_subsumed_hits_removed_overlaps_trimmed_dict.items():
                SeqIO.write(value['hit_sequence'], fasta_handle_nucl, 'fasta')

    def _create_stitched_contig(self):
        """
        Takes a dictionary of filtered, sorted, and trimmed hits, and creates a stitched contig SeqRecord by
        concatenating the corresponding sequences. If only one hit is present, return a SeqRecord of this hit.

        :return Bio.SeqRecord.SeqRecord: no_stitched_contig or stitched_contig, depending on number of hits
        """

        if not self.hits_subsumed_hits_removed_overlaps_trimmed_dict:
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
                                f'single BLASTn hit.'
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
                        f'BLASTn hit.'
            self._write_genes_with_stitched_contig(log_entry)

        return stitched_contig_seqrecord

    def _no_stitched_contig(self):
        """
        Identifies the single longest BLASTn hit; does not attempt to stitch multiple hits together into a
        stitched contig.

        :return Bio.SeqRecord.SeqRecord: SeqRecord for the single longest Exonerate hit.
        """

        if not self.hits_filtered_by_pct_similarity_dict:
            return None

        sample_name = os.path.split(self.prefix)[-1]
        gene_name = os.path.split(self.prefix)[-2]

        # Don't overwrite self.hits_subsumed_hits_removed_dict:
        blast_hits_subsumed_hits_removed_copy = copy.deepcopy(self.hits_subsumed_hits_removed_dict)

        sorted_by_hit_length = sorted(blast_hits_subsumed_hits_removed_copy.values(),
                                      key=lambda x: len(x['hit_sequence']), reverse=True)

        if len(sorted_by_hit_length) == 0:
            raise ValueError(f'The list sorted_by_hit_length for gene {gene_name} is empty!')

        sorted_by_hit_length[0]['hit_sequence'].description = f'Flag no_stitched_contig used. Single longest hit ' \
                                                              f'{sorted_by_hit_length[0]["hit_sequence"].id}'
        sorted_by_hit_length[0]['hit_sequence'].id = sample_name
        sorted_by_hit_length[0]['hit_sequence'].name = sample_name

        # Write report file:
        log_entry = f'{sample_name},{gene_name}, Stitched contig step skipped (user provided the ' \
                    f'--no_stitched_contig flag to hybpiper assemble). Gene sequence contains the longest BLASTn ' \
                    f'hit sequence only.'
        self._write_genes_with_stitched_contig(log_entry)

        return sorted_by_hit_length[0]['hit_sequence']

    def write_stitched_contig_to_file(self):
        """
        Writes DNA (suffix '.FNA') sequence for the stitched_contig (or single remaining hit sequence) to fasta file.

        :return NoneType: no explicit return
        """

        gene_name = os.path.split(self.prefix)[-2]
        dna_seqrecord_to_write = self.stitched_contig_seqrecord

        with open(f'{self.prefix}/sequences/FNA/{gene_name}.FNA', 'w') as fna_handle:
            SeqIO.write(dna_seqrecord_to_write, fna_handle, 'fasta')

    def _recover_long_paralogs(self):
        """
        Determines whether there are multiple BLASTn hits that are >75% the length of the query protein. If so,
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
        with an equal high similarity percentage, the annotated paralog will be the first in the list returned by
        Python's list.sort() method.

        :param collections.defaultdict paralog_dicts: dictionary of dictionaries; hitname: {key:value hit data}
        :return collections.defaultdict: paralog_dicts, with the highest similarity paralog annotated in seq.description
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

    def write_exonerate_stats_file(self):
        """
        Write a *.tsv file with stats from the initial BLASTn run following filtering by similarity,
        and stats for subsequent filtering/trimming steps.

        :return NoneType: no explicit return
        """

        headers = ['query_id',
                   'query_length',
                   'hit_id',
                   'query_HSP_range_limits_original',
                   'query_HSP_range_limits_trimmed',
                   'query_HSPFragment_ranges',
                   'hit_percent_similarity_original',
                   'hit_percent_similarity_trimmed',
                   'hit_strand',
                   'hit_HSP_range_limits_original',
                   'hit_HSP_range_limits_trimmed',
                   'hit_HSPFragment_ranges_original',
                   'hit_HSPFragment_ranges_trimmed',
                   '3-prime_bases_trimmed']

        if not self.hits_filtered_by_pct_similarity_dict:
            return

        def nested_dict_iterator(nested_dict):
            """
            Nested function to yield the key and value of a given nested dict.

            :param dict, collections.defaultdict nested_dict:
            :return: yields the key and value of the nested dict
            """

            for dict_key, dict_value in nested_dict.items():
                if isinstance(dict_value, dict):
                    yield dict_key, dict_value

        with open(f'{self.prefix}/blast_stats.tsv', 'w') as blast_stats_handle:
            headers = '\t'.join(headers)
            blast_stats_handle.write(f"\t{headers}\n")

            # Write stats for hits filtered by similarity:
            nested_dict_iter = nested_dict_iterator(self.hits_filtered_by_pct_similarity_dict)
            first_hit_dict_key, first_hit_dict_value = next(nested_dict_iter)  # get first line
            blast_stats_handle.write(f"Hits filtered > {self.similarity_threshold} percent similarity"
                                     f"\t{self.query_id}"
                                     f"\t{self.query_length}"
                                     f"\t{str(first_hit_dict_key.split(',')[0])}"
                                     f"\t{first_hit_dict_value['query_range_original']}"
                                     f"\t{first_hit_dict_value['query_range']}"
                                     f"\t{first_hit_dict_value['query_range_all']}"
                                     f"\t{first_hit_dict_value['hit_similarity_original']}"
                                     f"\t{first_hit_dict_value['hit_similarity']}"
                                     f"\t{first_hit_dict_value['hit_strand']}"
                                     f"\t{first_hit_dict_value['hit_range_original']}"
                                     f"\t{first_hit_dict_value['hit_range']}"
                                     f"\t{first_hit_dict_value['hit_range_all_original']}"
                                     f"\t{first_hit_dict_value['hit_range_all']}"
                                    f"\tN/A\n")  # as no trimming performed yet
            for key, value in nested_dict_iter:  # Write remaining lines
                blast_stats_handle.write(f"\t{self.query_id}"
                                         f"\t{self.query_length}"
                                         f"\t{str(key.split(',')[0])}"
                                         f"\t{value['query_range_original']}"
                                         f"\t{value['query_range']}"
                                         f"\t{value['query_range_all']}"
                                         f"\t{value['hit_similarity_original']}"
                                         f"\t{value['hit_similarity']}"
                                         f"\t{value['hit_strand']}"
                                         f"\t{value['hit_range_original']}"
                                         f"\t{value['hit_range']}"
                                         f"\t{value['hit_range_all_original']}"
                                         f"\t{value['hit_range_all']}"
                                         f"\tN/A\n")

            # Write stats for hits filtered by similarity and with subsumed hits removed:
            nested_dict_iter = nested_dict_iterator(self.hits_subsumed_hits_removed_dict)
            first_hit_dict_key, first_hit_dict_value = next(nested_dict_iter)  # get first line
            blast_stats_handle.write(f"Hits with subsumed hits removed"
                                     f"\t{self.query_id}"
                                     f"\t{self.query_length}"
                                     f"\t{str(first_hit_dict_key.split(',')[0])}"
                                     f"\t{first_hit_dict_value['query_range_original']}"
                                     f"\t{first_hit_dict_value['query_range']}"
                                     f"\t{first_hit_dict_value['query_range_all']}"
                                     f"\t{first_hit_dict_value['hit_similarity_original']}"
                                     f"\t{first_hit_dict_value['hit_similarity']}"
                                     f"\t{first_hit_dict_value['hit_strand']}"
                                     f"\t{first_hit_dict_value['hit_range_original']}"
                                     f"\t{first_hit_dict_value['hit_range']}"
                                     f"\t{first_hit_dict_value['hit_range_all_original']}"
                                     f"\t{first_hit_dict_value['hit_range_all']}"
                                     f"\tN/A\n")  # as no trimming performed yet
            for key, value in nested_dict_iter:  # Write remaining lines
                blast_stats_handle.write(f"\t{self.query_id}"
                                         f"\t{self.query_length}"
                                         f"\t{str(key.split(',')[0])}"
                                         f"\t{value['query_range_original']}"
                                         f"\t{value['query_range']}"
                                         f"\t{value['query_range_all']}"
                                         f"\t{value['hit_similarity_original']}"
                                         f"\t{value['hit_similarity']}"
                                         f"\t{value['hit_strand']}"
                                         f"\t{value['hit_range_original']}"
                                         f"\t{value['hit_range']}"
                                         f"\t{value['hit_range_all_original']}"
                                         f"\t{value['hit_range_all']}"
                                         f"\tN/A\n")

            # Write stats for hits filtered by similarity and subsumed hits removed and trimmed:
            nested_dict_iter = nested_dict_iterator(self.hits_subsumed_hits_removed_overlaps_trimmed_dict)
            first_hit_dict_key, first_hit_dict_value = next(nested_dict_iter)  # get first line
            trimmed_bases = first_hit_dict_value['hit_sequence'].description.split(':')[-1].strip()
            blast_stats_handle.write(f"Hits with subsumed hits removed and overlaps trimmed"
                                     f"\t{self.query_id}"
                                     f"\t{self.query_length}"
                                     f"\t{str(first_hit_dict_key.split(',')[0])}"
                                     f"\t{first_hit_dict_value['query_range_original']}"
                                     f"\t{first_hit_dict_value['query_range']}"
                                     f"\t{first_hit_dict_value['query_range_all']}"
                                     f"\t{first_hit_dict_value['hit_similarity_original']}"
                                     f"\t{first_hit_dict_value['hit_similarity']}"
                                     f"\t{first_hit_dict_value['hit_strand']}"
                                     f"\t{first_hit_dict_value['hit_range_original']}"
                                     f"\t{first_hit_dict_value['hit_range']}"
                                     f"\t{first_hit_dict_value['hit_range_all_original']}"
                                     f"\t{first_hit_dict_value['hit_range_all']}"
                                     f"\t{trimmed_bases}\n")
            for key, value in nested_dict_iter:  # Write remaining lines
                trimmed_bases = value['hit_sequence'].description.split()[-1]
                blast_stats_handle.write(f"\t{self.query_id}"
                                         f"\t{self.query_length}"
                                         f"\t{str(key.split(',')[0])}"
                                         f"\t{value['query_range_original']}"
                                         f"\t{value['query_range']}"
                                         f"\t{value['query_range_all']}"
                                         f"\t{value['hit_similarity_original']}"
                                         f"\t{value['hit_similarity']}"
                                         f"\t{value['hit_strand']}"
                                         f"\t{value['hit_range_original']}"
                                         f"\t{value['hit_range']}"
                                         f"\t{value['hit_range_all_original']}"
                                         f"\t{value['hit_range_all']}"
                                         f"\t{trimmed_bases}\n")

            # Write stats for paralogs:
            if self.long_paralogs_dict:
                dict_iter = iter((key, value) for key, value in self.long_paralogs_dict.items())
                first_hit_dict_key, first_hit_dict_value = next(dict_iter)  # get first line
                blast_stats_handle.write(
                    (f"Hits corresponding to long paralogs"
                     f"\t{self.query_id}"
                     f"\t{self.query_length}"
                     f"\t{str(first_hit_dict_key.split(',')[0])}"
                     f"\t{first_hit_dict_value['query_range_original']}"
                     f"\t{first_hit_dict_value['query_range']}"
                     f"\t{first_hit_dict_value['query_range_all']}"
                     f"\t{first_hit_dict_value['hit_similarity_original']}"
                     f"\t{first_hit_dict_value['hit_similarity']}"
                     f"\t{first_hit_dict_value['hit_strand']}"
                     f"\t{first_hit_dict_value['hit_range_original']}"
                     f"\t{first_hit_dict_value['hit_range']}"
                     f"\t{first_hit_dict_value['hit_range_all_original']}"
                     f"\t{first_hit_dict_value['hit_range_all']}"
                     f"\tN/A\n"))
                for key, value in dict_iter:  # Write remaining lines
                    blast_stats_handle.write((f"\t{self.query_id}"
                                                  f"\t{self.query_length}"
                                                  f"\t{str(key.split(',')[0])}"
                                                  f"\t{value['query_range_original']}"
                                                  f"\t{value['query_range']}"
                                                  f"\t{value['query_range_all']}"
                                                  f"\t{value['hit_similarity_original']}"
                                                  f"\t{value['hit_similarity']}"
                                                  f"\t{value['hit_strand']}"
                                                  f"\t{value['hit_range_original']}"
                                                  f"\t{value['hit_range']}"
                                                  f"\t{value['hit_range_all_original']}"
                                                  f"\t{value['hit_range_all']}"
                                                  f"\tN/A\n"))

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

        if not self.hits_subsumed_hits_removed_overlaps_trimmed_dict:
            return None

        hit_ranges_dict = defaultdict(list)

        cumulative_hit_length = 0  # track to adjust coordinates of hits to match exon-only stitched contig sequence
        for hit, hit_dict_values in self.hits_subsumed_hits_removed_overlaps_trimmed_dict.items():
            if self.verbose_logging:
                self.logger.debug(f'hit is: {hit}')
            spades_name = hit.split(',')[0]
            raw_spades_contig_length = len(self.spades_assembly_dict[spades_name])
            hit_exonerate_sequence_length = len(hit_dict_values['hit_sequence'].seq)
            hit_range_all = hit_dict_values['hit_range_all']
            hit_inter_ranges = hit_dict_values['hit_inter_ranges']

            if len(hit_inter_ranges) == 0:  # i.e. no introns
                if self.verbose_logging:
                    self.logger.debug(f'len(hit_inter_ranges) for hit {hit} is 0, no intron coordinates for chimera test')
                hit_ranges_dict[hit].append('no introns')
                cumulative_hit_length += hit_exonerate_sequence_length
            else:
                if hit_dict_values['hit_strand'] == -1:  # Convert ranges so that they apply to the revcomp contig
                    if self.verbose_logging:
                        self.logger.debug(f'hit_inter_ranges before conversion is: {hit_inter_ranges}')
                    hit_inter_ranges = self.convert_coords_revcomp(hit_inter_ranges, raw_spades_contig_length)
                    if self.verbose_logging:
                        self.logger.debug(f'hit_inter_ranges after conversion is: {hit_inter_ranges}')
                        self.logger.debug(f'hit_range_all before conversion is: {hit_range_all}')
                    hit_range_all = self.convert_coords_revcomp(hit_range_all, raw_spades_contig_length)
                    if self.verbose_logging:
                        self.logger.debug(f'hit_range_all after conversion is: {hit_range_all}')
                else:
                    if self.verbose_logging:
                        self.logger.debug(f'hit_inter_ranges is: {hit_inter_ranges}')
                        self.logger.debug(f'hit_range_all is: {hit_range_all}')

                # Adjust range coordinates so that they start at zero (i.e. first nucleotide position of the Exonerate
                # fasta sequence for this hit:
                hit_start_coordinate = hit_range_all[0][0]
                hit_range_all_start_at_zero = [(item[0] - hit_start_coordinate, item[1] - hit_start_coordinate)
                                               for item in hit_range_all]
                if self.verbose_logging:
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

        if self.verbose_logging:
            self.logger.debug(f'hit_ranges_dict is: {hit_ranges_dict}')

        return hit_ranges_dict

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

        if self.verbose_logging:
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
                        f'-Xmx{self.chimera_bbmap_memory}m ' \
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

        except subprocess.CalledProcessError as exc:
            fill = textwrap.fill(f'{"[ERROR]:":10} Running bbmap.sh during the stitched-contig chimera test failed '
                                 f'for gene {gene_name}. No stitched-contig chimera test will be performed. See the '
                                 f'*.log file in the sample directory for error details.', width=90,
                                 subsequent_indent=' ' * 11)
            self.logger.info(f'\n{fill}')
            self.logger.debug(f'bbmap_command FAILED. Output is: {exc}')
            self.logger.debug(f'bbmap_command stdout is: {exc.stdout}')
            self.logger.debug(f'bbmap_command stderr is: {exc.stderr}')
            return False

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

        if self.verbose_logging:
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

        if self.verbose_logging:
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

            if self.verbose_logging:
                self.logger.debug(f'forward_enclosing_contig_range is: {forward_enclosing_contig_range}')
                self.logger.debug(f'reverse_enclosing_contig_range is: {reverse_enclosing_contig_range}')

            # Make sure each read is assigned to a contig range, or skip if it overlaps contig join:
            if not forward_enclosing_contig_range or not reverse_enclosing_contig_range:
                if self.verbose_logging:
                    self.logger.debug(f'One or both read pairs do not map within the range of a single contig (i.e. '
                                      f'they overlap a coordinate where two contigs have been concatenated). Skipping '
                                      f'read pair')
                continue

            forward_contig_exon_ranges = self.stitched_contig_hit_ranges[forward_enclosing_contig_name]
            reverse_contig_exon_ranges = self.stitched_contig_hit_ranges[reverse_enclosing_contig_name]

            self.logger.debug(f'forward_contig_exon_ranges is: {forward_contig_exon_ranges}')
            self.logger.debug(f'reverse_contig_exon_ranges is: {reverse_contig_exon_ranges}')

            # Check if the SPAdes only name each contig are the same:
            if forward_enclosing_contig_name_spades_only == reverse_enclosing_contig_name_spades_only:
                if self.verbose_logging:
                    self.logger.debug(f'Both reads occur in the same SPAdes contig. Skipping read pair!')
                continue
            else:  # If contig is not the same, check that reads don't overlap exon-intron boundaries within each contig
                if self.verbose_logging:
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
                    if self.verbose_logging:
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

    def write_no_stitched_contig(self):
        """
        Writes DNA (suffix '.FNA') sequence for the single longest Exonerate hit (or single remaining hit sequence)
        to fasta file.

        :return NoneType: no explicit return
        """

        gene_name = os.path.split(self.prefix)[-2]
        hit_id = self.no_stitched_contig_seqrecord.id
        name = self.no_stitched_contig_seqrecord.name
        description = self.no_stitched_contig_seqrecord.description

        dna_seqrecord_to_write = self.no_stitched_contig_seqrecord

        with open(f'{self.prefix}/sequences/FNA/{gene_name}.FNA', 'w') as fna_handle:
            SeqIO.write(dna_seqrecord_to_write, fna_handle, 'fasta')

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
        interleaved_reads = f'{gene_folder}/{gene_folder}_interleaved.fasta'  # FIXME same for single end?

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


def parse_spades_and_best_reference(assemblyfile, locus_file, prefix):
    """
    Return a SeqIO dictionary for the SPAdes contigs (assemblyfile) and the 'best' locus
    reference for the sample/gene.

    :param str assemblyfile: name of the FASTA file containing DNA sequence assembly
    :param str locus_file: name of the FASTA file containing the target locus sequence
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
        locus_file = open(locus_file)
    except IOError:
        logger.debug(f'The file {locus_file} could not be opened!')
        return

    spades_assembly_dict = SeqIO.to_dict(SeqIO.parse(assemblyfile, 'fasta'))
    best_locus_ref_dict = SeqIO.to_dict(SeqIO.parse(locus_file, 'fasta'))

    assemblyfile.close()
    locus_file.close()

    return spades_assembly_dict, best_locus_ref_dict


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
    parser.add_argument("--bbmap_memory", help="MB memory (RAM) to use for bbmap.sh", default=250, type=int)
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
    parser.add_argument("--no_pad_stitched_contig_gaps_with_n",
                        help="When constructing stitched contigs, do not pad any gaps between hits (with respect to "
                             "the 'best' protein reference) with a number of Ns corresponding to the "
                             "reference gap multiplied by 3. Default is %(default)s.",
                        action="store_false",
                        dest='stitched_contig_pad_n',
                        default=True)
    parser.add_argument('--chimeric_stitched_contig_check',
                        help='Attempt to determine whether a stitched contig is a potential chimera of contigs from '
                             'multiple paralogs. Default is %(default)s.',
                        action='store_true',
                        dest='chimera_check',
                        default=False)
    parser.add_argument("--no_intronerate",
                        help="Do no run intronerate to recover fasta files for supercontig with introns (if present), "
                             "and introns-only", action="store_true", dest='no_intronerate', default=False)
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
    parser.add_argument('--exonerate_hit_sliding_window_size',
                        help='Size of the sliding window (in amino-acids) when trimming termini of Exonerate '
                             'hits. Default is %(default)s.',
                        default=3,
                        type=int)
    parser.add_argument('--exonerate_hit_sliding_window_thresh',
                        help='Percentage similarity threshold for the sliding window (in amino-acids) when trimming '
                             'termini of Exonerate hits. Default is %(default)s.',
                        default=55,
                        type=int)
    parser.add_argument('--exonerate_skip_hits_with_frameshifts',
                        help='Skip Exonerate hits where the SPAdes sequence contains a frameshift. Default is '
                             '%(default)s.',
                        action='store_true',
                        dest='skip_frameshifts',
                        default=False)
    parser.add_argument('--exonerate_skip_hits_with_internal_stop_codons',
                        help='Skip Exonerate hits where the SPAdes sequence contains an internal in-frame stop codon. '
                             'See:\nhttps://github.com/mossmatters/HybPiper/wiki/Troubleshooting,-common-issues,'
                             '-and-recommendations#31-sequences-containing-stop-codons.\n A single terminal stop codon '
                             'is allowed, but see option "--exonerate_skip_hits_with_terminal_stop_codons" below. '
                             'Default is %(default)s.',
                        action='store_true',
                        dest='skip_internal_stops',
                        default=False)
    parser.add_argument('--exonerate_skip_hits_with_terminal_stop_codons',
                        help='Skip Exonerate hits where the SPAdes sequence contains a single terminal stop codon. '
                             'Only applies when option "--exonerate_skip_hits_with_internal_stop_codons" is also '
                             'provided. Only use this flag if your target file exclusively contains protein-coding '
                             'genes with no stop codons included, and you would like to prevent any in-frame stop '
                             'codons in the output sequences. Default is %(default)s.',
                        action='store_true',
                        dest='skip_terminal_stops',
                        default=False)
    parser.add_argument('--verbose_logging',
                        help='If supplied, enable verbose login. NOTE: this can increase the size of the log '
                             'files by an order of magnitude.',
                        action='store_true',
                        dest='verbose_logging',
                        default=False)

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

    exonerate_result = parse_exonerate_and_get_stitched_contig(
        exonerate_text_output,
        query_file=args.proteinfile,
        paralog_warning_min_length_percentage=
        args.paralog_warning_min_length_percentage,
        thresh=args.thresh,
        logger=logger,
        prefix=prefix,
        chimera_check=args.chimera_check,
        discordant_cutoff=
        args.chimeric_stitched_contig_discordant_reads_cutoff,
        edit_distance=args.chimeric_stitched_contig_edit_distance,
        bbmap_subfilter=args.bbmap_subfilter,
        bbmap_memory=args.bbmap_memory,
        bbmap_threads=args.bbmap_threads,
        interleaved_fasta_file=path_to_interleaved_fasta,
        no_stitched_contig=args.no_stitched_contig,
        stitched_contig_pad_n=args.stitched_contig_pad_n,
        spades_assembly_dict=spades_assembly_dict,
        depth_multiplier=args.depth_multiplier,
        keep_intermediate_files=args.keep_intermediate_files,
        exonerate_hit_sliding_window_size=args.exonerate_hit_sliding_window_size,
        exonerate_hit_sliding_window_thresh=args.exonerate_hit_sliding_window_thresh,
        exonerate_skip_frameshifts=args.skip_frameshifts,
        exonerate_skip_internal_stops=args.skip_internal_stops,
        exonerate_skip_terminal_stops=args.skip_terminal_stops,
        verbose_logging=args.verbose_logging)

    if not exonerate_result.stitched_contig_seqrecord:
        return

    logger.debug(f'There were {len(exonerate_result.hits_filtered_by_pct_similarity_dict)} Exonerate '
                 f'hits for {args.proteinfile} after filtering by similarity threshold {args.thresh}.')


########################################################################################################################
# Run the script
########################################################################################################################
if __name__ == "__main__":
    standalone()

################################################## END OF SCRIPT #######################################################
