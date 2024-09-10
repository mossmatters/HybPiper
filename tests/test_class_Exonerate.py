"""
Tests for class Exonerate in module exonerate_hits.py.
"""

import os
import sys
import logging

try:
    import pyfakefs
except ImportError:
    sys.exit(f'The required Python package pyfakefs (https://github.com/pytest-dev/pyfakefs) is not found.\n'
             f'Is it installed for the Python installation used to run this test?')

from pyfakefs.fake_filesystem_unittest import TestCase
from Bio.SearchIO import parse
from hybpiper.exonerate_hits import Exonerate
from hybpiper.exonerate_hits import parse_spades_and_best_reference
from hybpiper.exonerate_hits import create_output_directories

TEST_DIR_EXONERATE = "Exonerate"
TEST_DIR_SPADES = "SPAdes_assemblies"
TEST_DIR_TARGET = "Targets"
TEST_DIR_READS = "Read_files"


class ExonerateClass(TestCase):

    def setUp(self):
        self.setUpPyfakefs()
        self.fs.add_real_directory(TEST_DIR_EXONERATE)
        self.fs.add_real_directory(TEST_DIR_SPADES)
        self.fs.add_real_directory(TEST_DIR_TARGET)
        self.fs.add_real_directory(TEST_DIR_READS)

    def test_processing_exonerate_results_frameshifts_qresult(self):
        """
        Test processing of an Exonerate qresult object parsed with Bio.SearchIO.parse (format = exonerate-text) and
        generation of the correct stitched-contig *.FNA sequence.
        """

        # Create a logger:
        logger = logging.getLogger()

        # Set directory names, and prefix for instance of class Exonerate:
        sample_name = 'cocos_nucifera_006118'
        gene_name = '4527'
        prefix_exonerate = '4527/cocos_nucifera_006118'

        # Create required directories (usually done by assemble.py):
        os.makedirs(f'{sample_name}')
        os.chdir(sample_name)  # move to sample directory
        os.mkdir(gene_name)

        # Get the target file protein query fasta file:
        proteinfile = os.path.join('..', TEST_DIR_TARGET, '4527_target.fasta')

        # Get the SPAdes contigs fasta file:
        assemblyfile = os.path.join('..', TEST_DIR_SPADES, '4527_contigs.fasta')

        # Get the interleaved read file:
        interleaved_fasta_file = os.path.join('..', TEST_DIR_READS, '4527_interleaved.fasta')

        # Get the Exonerate output text file:
        exn_file = os.path.join('..', TEST_DIR_EXONERATE, 'exonerate_results_frameshifts.fasta')

        # Create sequence output directories:
        create_output_directories(prefix_exonerate, assemblyfile)

        # Parse the Exonerate text output file with SearchIO:)
        exonerate_hits_from_alignment = list(parse(exn_file, 'exonerate-text'))

        # Read the SPAdes contigs and the 'best' protein reference seq into SeqIO dictionaries:
        spades_assembly_dict, best_protein_ref_dict = parse_spades_and_best_reference(assemblyfile,
                                                                                      proteinfile,
                                                                                      prefix_exonerate)

        self.assertEqual(
            "[QueryResult(id='Arecales_Arecaceae_Chamaedorea_seifrizii_Timilsena_et_al-4527', 5 hits)]",
            str(exonerate_hits_from_alignment)
        )

        # Create Exonerate class instance - writes to pyfakefs:
        exonerate_result = Exonerate(searchio_object=exonerate_hits_from_alignment,
                                     query_file=proteinfile,
                                     paralog_warning_min_length_percentage=0.75,
                                     thresh=55,
                                     chimera_check=False,
                                     logger=logger,
                                     prefix=prefix_exonerate,
                                     discordant_cutoff=100,
                                     edit_distance=7,
                                     bbmap_subfilter=7,
                                     bbmap_memory=250,
                                     bbmap_threads=2,
                                     interleaved_fasta_file=interleaved_fasta_file,
                                     no_stitched_contig=False,
                                     spades_assembly_dict=spades_assembly_dict,
                                     depth_multiplier=10,
                                     keep_intermediate_files=False,
                                     trim_hit_sliding_window_size=5,
                                     trim_hit_sliding_window_thresh=75,
                                     exonerate_skip_frameshifts=False,
                                     verbose_logging=False)

        exonerate_result.write_stitched_contig_to_file()  # Writes to pyfakefs

        with open(f'{prefix_exonerate}/sequences/FNA/4527.FNA') as fna_handle:
            fasta_seq = fna_handle.read()

        self.assertEqual('>cocos_nucifera_006118 multi_hit_stitched_contig_comprising_3_hits\n'
                         'GAGAGGGTAGTGGTGCTGGTCATTGGTGGTGGAGGAAGGGAGCATGCACTTTGCTATGCC\n'
                         'TTGAAGCAGTCTCCATCCTGTGATGCAGTTTTTTGTGCCCCTGGTAATGCAGGGATTGGT\n'
                         'CAATCTGGGGATGCCACCTGCATATCAGACCTAAATATCTCCGACAGTGCAGCTGTGATA\n'
                         'TCCTTTTGCCGCATGTGGGGAGTTGGCCTGGTGGTAATTGGTCCAGAAGCTCCCTTAGTT\n'
                         'GCTGGCCTTGCGAATGACCTCGTTAAGGCTGGGATCCCAACGTTTGGTCCATCAGCAGAG\n'
                         'GCTGCAGCATTGGAGGGATCGAAGGACTTCATGAAGAAGTTGTGCGACAAATATGGCATC\n'
                         'CCTACAGCACAGTATCAAACCTTTACAAACCCTTCTGATGCAAAGCAATACGTCAAGGAG\n'
                         'CAAGGAGCACCTATTGTTGTGAAGGCTGACGGATTGGCTGCTGGGAAGGGAGTCATTGTT\n'
                         'GCCATGACTTTGGAGGAGGCATATGAAGCTATAGACTCGATGCTGGTTAATGGTGCTTTT\n'
                         'GGATCTGCTGGTTCCCGAGTCATTATTGAAGAATTCCTAGAGGGTGAAGAAGCTTCTTTC\n'
                         'TTTGCTCTCGTCGATGGTGAAACTGCCTTACCGCTTGAATCTGCACAGGACCATAAAAGA\n'
                         'GTGGGTGATGGAGATACTGGCCCAAATACAGGTGGAATGGGTGCATATTCTCCTGCACCG\n'
                         'GTACTAACAAAAGAACTTCAGTCTGTTGTCATGGAATCAATAATCCTCCCAACAGTGAAG\n'
                         'GGGATGGCAGCTGAAGGGTGCAAATTTGTTGGTGTGTTGTATGCTGGGCTTATGATTGAG\n'
                         'AAAAAATCTAGGCTACCTAAGCTCATTGAGTACAATGTCCGCTTTGGAGATCCTGAATGC\n'
                         'CAGGTTCTGATGATGCGGTTGGAGTCCGATCTGGCACAGGTTCTGCTTGCAGCTTCTCGG\n'
                         'GGAGAGTTGGGCAATGCGTCGCTGAGTTGGTCACCTGGACCAGCTATGGTGGTTGTGATG\n'
                         'GCTAGTCGAGGTTATCCAGGTGCCTATGAGAAGGGCACTGTGATAAGAAACCTAGAAGAA\n'
                         'GCAGAGCTGATTTCTCCGATGGTTAAGATATTTCATGCTGGAACAGCTCTGGACTCAAAT\n'
                         'GGCAATTTCATTGCCTCTGGGGGTCGTGTGCTTGGGGTAACCGCAAAGGGGAAAGACATA\n'
                         'GAAGAAGCGAGAACAAGGGCTTATGATGCAGTTGAAGCAATCGACTGGCCTGAAGGATTC\n'
                         'TATAGGCATGATATTGGCTGGAGAGCA\n',
                         fasta_seq)

    def test_processing_exonerate_results_frameshifts_qresult_skip_frameshifts(self):
        """
        Test processing of an Exonerate qresult object parsed with Bio.SearchIO.parse (format = exonerate-text),
        using option exonerate_skip_frameshifts = True, and generation of the correct stitched-contig *.FNA sequence.
        """

        # Create a logger:
        logger = logging.getLogger()

        # Set directory names, and prefix for instance of class Exonerate:
        sample_name = 'cocos_nucifera_006118'
        gene_name = '4527'
        prefix_exonerate = '4527/cocos_nucifera_006118'

        # Create required directories (usually done by assemble.py):
        os.makedirs(f'{sample_name}')
        os.chdir(sample_name)  # move to sample directory
        os.mkdir(gene_name)

        # Get the target file protein query fasta file:
        proteinfile = os.path.join('..', TEST_DIR_TARGET, '4527_target.fasta')

        # Get the SPAdes contigs fasta file:
        assemblyfile = os.path.join('..', TEST_DIR_SPADES, '4527_contigs.fasta')

        # Get the interleaved read file:
        interleaved_fasta_file = os.path.join('..', TEST_DIR_READS, '4527_interleaved.fasta')

        # Get the Exonerate output text file:
        exn_file = os.path.join('..', TEST_DIR_EXONERATE, 'exonerate_results_frameshifts.fasta')

        # Create sequence output directories:
        create_output_directories(prefix_exonerate, assemblyfile)

        # Parse the Exonerate text output file with SearchIO:)
        exonerate_hits_from_alignment = list(parse(exn_file, 'exonerate-text'))

        # Read the SPAdes contigs and the 'best' protein reference seq into SeqIO dictionaries:
        spades_assembly_dict, best_protein_ref_dict = parse_spades_and_best_reference(assemblyfile,
                                                                                      proteinfile,
                                                                                      prefix_exonerate)

        self.assertEqual(
            "[QueryResult(id='Arecales_Arecaceae_Chamaedorea_seifrizii_Timilsena_et_al-4527', 5 hits)]",
            str(exonerate_hits_from_alignment)
        )

        # Create Exonerate class instance - writes to pyfakefs:
        exonerate_result = Exonerate(searchio_object=exonerate_hits_from_alignment,
                                     query_file=proteinfile,
                                     paralog_warning_min_length_percentage=0.75,
                                     thresh=55,
                                     chimera_check=False,
                                     logger=logger,
                                     prefix=prefix_exonerate,
                                     discordant_cutoff=100,
                                     edit_distance=7,
                                     bbmap_subfilter=7,
                                     bbmap_memory=250,
                                     bbmap_threads=2,
                                     interleaved_fasta_file=interleaved_fasta_file,
                                     no_stitched_contig=False,
                                     spades_assembly_dict=spades_assembly_dict,
                                     depth_multiplier=10,
                                     keep_intermediate_files=False,
                                     trim_hit_sliding_window_size=5,
                                     trim_hit_sliding_window_thresh=75,
                                     exonerate_skip_frameshifts=True,
                                     verbose_logging=False)

        exonerate_result.write_stitched_contig_to_file()  # Writes to pyfakefs

        with open(f'{prefix_exonerate}/sequences/FNA/4527.FNA') as fna_handle:
            fasta_seq = fna_handle.read()

        self.assertEqual('>cocos_nucifera_006118 multi_hit_stitched_contig_comprising_2_hits\n'
                         'GAGAGGGTAGTGGTGCTGGTCATTGGTGGTGGAGGAAGGGAGCATGCACTTTGCTATGCC\n'
                         'TTGAAGCAGTCTCCATCCTGTGATGCAGTTTTTTGTGCCCCTGGTAATGCAGGGATTGGT\n'
                         'CAATCTGGGGATGCCACCTGCATATCAGACCTAAATATCTCCGACAGTGCAGCTGTGATA\n'
                         'TCCTTTTGCCGCATGTGGGGAGTTGGCCTGGTGGTAATTGGTCCAGAAGCTCCCTTAGTT\n'
                         'GCTGGCCTTGCGAATGACCTCGTTAAGGCTGGGATCCCAACGTTTGGTCCATCAGCAGAG\n'
                         'GCTGCAGCATTGGAGGGATCGAAGGACTTCATGAAGAAGTTGTGCGACAAATATGGCATC\n'
                         'CCTACAGCAAAGTATCAAACCTTTACAAACCCTTCTGATGCAAAGCAATACGTCAAGGAG\n'
                         'CAAGGAGCACCTATTGTTGTGAAGGCTGACGGATTGGCTGCTGGGAAGGGAGTCATTGTT\n'
                         'GCCATGACTTTGGAGGAGGCATATGAAGCTATAGACTCGATGCTGGTTAATGGTGCTTTT\n'
                         'GGATCTGCTGGTTCCCGAGTCATTATTGAAGAATTCCTAGAGGGTGAAGAAGCTTCTTTC\n'
                         'TTTGCTCTCGTCGATGGTGAAACTGCCTTACCGCTTGAATCTGCACAGGACCATAAAAGA\n'
                         'GTGGGTGATGGAGATACTGGCCCAAATACAGGTGGAATGGGTGCATATTCTCCTGCACCG\n'
                         'GTACTAACAAAAGAACTTCAGTCTGTTGTCATGGAATCAATAATCCTCCCAACAGTGAAG\n'
                         'GGGATGGCAGCTGAAGGGTGCAAATTTGTTGGTGTGTTGTATGCTGGGCTTATGATTGAG\n'
                         'AAAAAATCTAGGCTACCTAAGCTCATTGAGTACAATGTCCGCTTTGGAGATCCTGAATGC\n'
                         'CAGGTTCTGATGATGCGGTTGGAGTCCGATCTGGCACAGGTTCTGCTTGCAGCTTCTCGG\n'
                         'GGAGAGTTGGGCAATGCGTCGCTGAGTTGGTCACCTGGACCAGCTATGGTGGTTGTGATG\n'
                         'GCTAGTCGAGGTTATCCAGGTGCCTATGAGAAGGGCACTGTGATAAGAAACCTAGAAGAA\n'
                         'GCAGAGCTGATTTCTCCGATGGTTAAGATATTTCATGCTGGAACAGCTCTGGACTCAAAT\n'
                         'GGCAATTTCATTGCCTCTGGGGGTCGTGTGCTTGGGGTAACCGCAAAGGGGAAAGACATA\n'
                         'GAAGAAGCGAGAACAAGGGCTTATGATGCAGTTGAAGCAATCGACTGGCCTGAAGGATTC\n'
                         'TATAGGCATGATATTGGCTGGAGAGCA\n',
                         fasta_seq)
