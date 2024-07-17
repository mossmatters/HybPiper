#!/usr/bin/env python

"""
Contains subparsers used in hybpiper_main.py
"""

import logging
from hybpiper import utils

# Create logger:
logger = logging.getLogger(f'hybpiper.hybpiper_main.{__name__}')


def add_assemble_parser(subparsers):
    """
    Parser for the main assembly stage of HybPiper i.e. assemble.py.

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_assemble = subparsers.add_parser('assemble',
                                            help='Assemble gene, intron, and supercontig sequences for a sample')
    parser_assemble.add_argument('--readfiles', '-r',
                                 nargs='+',
                                 required=True,
                                 help='One or more read files to start the pipeline. If exactly two are specified, '
                                      'will assume it is paired Illumina reads.')
    group_1 = parser_assemble.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna',
                         dest='targetfile_dna',
                         default=False,
                         help='FASTA file containing DNA target sequences for each gene. The fasta headers must '
                              'follow the naming convention: >TaxonID-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa',
                         dest='targetfile_aa',
                         default=False,
                         help='FASTA file containing amino-acid target sequences for each gene. The fasta headers '
                              'must follow the naming convention: >TaxonID-geneName')
    group_2 = parser_assemble.add_mutually_exclusive_group()
    group_2.add_argument('--bwa',
                         dest='bwa',
                         action='store_true',
                         default=False,
                         help='Use BWA to search reads for hits to target. Requires BWA and a target file that is '
                              'nucleotides! Default is: %(default)s')
    group_2.add_argument('--diamond',
                         dest='diamond',
                         action='store_true',
                         default=False,
                         help='Use DIAMOND instead of BLASTx. Default is: %(default)s')
    parser_assemble.add_argument('--diamond_sensitivity',
                                 choices=['mid-sensitive', 'sensitive', 'more-sensitive', 'very-sensitive',
                                          'ultra-sensitive'],
                                 default=False,
                                 help='Use the provided sensitivity for DIAMOND searches.')
    parser_assemble.add_argument('--start_from',
                                 choices=['map_reads', 'distribute_reads', 'assemble_reads', 'exonerate_contigs'],
                                 dest='start_from',
                                 default='map_reads',
                                 help='Start the pipeline from the given step. Note that this relies on the presence '
                                      'of output files for previous steps, produced by a previous run attempt. '
                                      'Default is: %(default)s')
    parser_assemble.add_argument('--end_with',
                                 choices=['map_reads', 'distribute_reads', 'assemble_reads', 'exonerate_contigs'],
                                 dest='end_with',
                                 default='exonerate_contigs',
                                 help='End the pipeline at the given step. Default is: %(default)s')
    parser_assemble.add_argument('--force_overwrite',
                                 action='store_true',
                                 default=False,
                                 help='Overwrite any output from a previous run for pipeline steps >= --start_from '
                                      'and <= --end_with. Default is: %(default)s')
    parser_assemble.add_argument('--cpu',
                                 type=int,
                                 help='Limit the number of CPUs. Default is to use all cores available minus one.')
    parser_assemble.add_argument('--skip_targetfile_checks',
                                 action='store_true',
                                 default=False,
                                 help='Skip the target file checks. Can be used if you are confident that your target '
                                      'file has no issues (e.g. if you have previously run "hybpiper '
                                      'check_targetfile". Default is: %(default)s'),
    parser_assemble.add_argument('--distribute_low_mem',
                                 action='store_true',
                                 default=False,
                                 help='Distributing and writing reads to individual gene directories will be 40-50 '
                                      'percent slower, but can use less memory/RAM with large input files (see wiki). '
                                      'Default is: %(default)s')
    parser_assemble.add_argument('--evalue',
                                 type=float,
                                 default=1e-4,
                                 help='e-value threshold for blastx hits. Default is: %(default)s')
    parser_assemble.add_argument('--max_target_seqs',
                                 type=int,
                                 default=10,
                                 help='Max target seqs to save in BLASTx search. Default is: %(default)s')
    parser_assemble.add_argument('--cov_cutoff',
                                 default=8,
                                 help='Coverage cutoff for SPAdes. Default is: %(default)s')
    parser_assemble.add_argument('--single_cell_assembly',
                                 action='store_true',
                                 dest='spades_single_cell',
                                 default=False,
                                 help='Run SPAdes assemblies in MDA (single-cell) mode. Default is: %(default)s')
    parser_assemble.add_argument('--kvals', nargs='+',
                                 default=None,
                                 help='Values of k for SPAdes assemblies. SPAdes needs to be compiled to handle '
                                      'larger k-values! Default is auto-detection by SPAdes.')
    parser_assemble.add_argument('--thresh',
                                 type=int,
                                 default=55,
                                 help='Percent identity threshold for retaining Exonerate hits. Default is 55, '
                                      'but increase this if you are worried about contaminant sequences.')
    parser_assemble.add_argument('--paralog_min_length_percentage',
                                 default=0.75,
                                 type=float,
                                 help='Minimum length percentage of a contig Exonerate hit vs reference protein '
                                      'length for a paralog warning and sequence to be generated. Default is: '
                                      '%(default)s')
    parser_assemble.add_argument('--depth_multiplier',
                                 default=10,
                                 type=int,
                                 help='Assign a long paralog as the "main" sequence if it has a coverage depth '
                                      '<depth_multiplier> times all other long paralogs. Set to zero to not use '
                                      'depth. Default is: %(default)s')
    parser_assemble.add_argument('--prefix',
                                 default=None,
                                 help='Directory name for pipeline output, default is to use the FASTQ file name.')
    parser_assemble.add_argument('--timeout_assemble',
                                 default=0,
                                 type=int,
                                 help='Kill long-running gene assemblies if they take longer than X percent of '
                                      'average.')
    parser_assemble.add_argument('--timeout_exonerate_contigs',
                                 default=120,
                                 type=int,
                                 help='Kill long-running processes if they take longer than X seconds. Default is: '
                                      '%(default)s')
    parser_assemble.add_argument('--target',
                                 default=None,
                                 help='Use the target file sequence with this taxon name in Exonerate searches for '
                                      'each gene. Other targets for that gene will be used only for read sorting. Can '
                                      'be a tab-delimited file (one <gene>\\t<taxon_name> per line) or a single taxon '
                                      'name.')
    parser_assemble.add_argument('--exclude',
                                 default=None,
                                 help='Do not use any sequence with the specified taxon name string in Exonerate '
                                      'searches. Sequenced from this taxon will still be used for read sorting.')
    parser_assemble.add_argument('--unpaired',
                                 default=False,
                                 help='Include a single FASTQ file with unpaired reads along with two paired read '
                                      'files')
    parser_assemble.add_argument('--no_stitched_contig',
                                 dest='no_stitched_contig',
                                 action='store_true',
                                 default=False,
                                 help='Do not create any stitched contigs. The longest single Exonerate hit will be '
                                      'used. Default is: %(default)s.')
    parser_assemble.add_argument('--no_pad_stitched_contig_gaps_with_n',
                                 action="store_false",
                                 dest='stitched_contig_pad_n',
                                 default=True,
                                 help='When constructing stitched contigs, do not pad any gaps between hits (with '
                                      'respect to the "best" protein reference) with a number of Ns corresponding to '
                                      'the reference gap multiplied by 3. Default is: %(default)s.')
    parser_assemble.add_argument('--chimeric_stitched_contig_check',
                                 action='store_true',
                                 dest='chimera_check',
                                 default=False,
                                 help='Attempt to determine whether a stitched contig is a potential '
                                      'chimera of contigs from multiple paralogs. Default is: %(default)s.')
    parser_assemble.add_argument('--bbmap_memory',
                                 default=1000,
                                 type=int,
                                 help='MB memory (RAM) to use for bbmap.sh with exonerate_hits.py. Default: is '
                                      '%(default)s.')
    parser_assemble.add_argument('--bbmap_subfilter',
                                 default=7,
                                 type=int,
                                 help='Ban alignments with more than this many substitutions. Default is: %(default)s.')
    parser_assemble.add_argument('--bbmap_threads',
                                 default=1,
                                 type=int,
                                 help='Number of threads to use for BBmap when searching for chimeric stitched contig. '
                                      'Default is: %(default)s.')
    parser_assemble.add_argument('--chimeric_stitched_contig_edit_distance',
                                 default=5,
                                 type=int,
                                 help='Minimum number of differences between one read of a read pair vs the '
                                      'stitched contig reference for a read pair to be flagged as discordant. Default '
                                      'is: %(default)s.')
    parser_assemble.add_argument('--chimeric_stitched_contig_discordant_reads_cutoff',
                                 default=5,
                                 type=int,
                                 help='Minimum number of discordant reads pairs required to flag a stitched contig as '
                                      'a potential chimera of contigs from multiple paralogs. Default is %(default)s.')
    parser_assemble.add_argument('--exonerate_hit_sliding_window_size',
                                 default=3,
                                 type=int,
                                 help='Size of the sliding window (in amino-acids) when trimming termini of Exonerate '
                                      'hits. Default is: %(default)s.')
    parser_assemble.add_argument('--exonerate_hit_sliding_window_thresh',
                                 default=55,
                                 type=int,
                                 help='Percentage similarity threshold for the sliding window (in amino-acids) when '
                                      'trimming termini of Exonerate hits. Default is: %(default)s.')
    parser_assemble.add_argument('--exonerate_skip_hits_with_frameshifts',
                                 action='store_true',
                                 dest='skip_frameshifts',
                                 default=False,
                                 help='Skip Exonerate hits where the SPAdes sequence contains a frameshift. '
                                      'See:\nhttps://github.com/mossmatters/HybPiper/wiki/Troubleshooting-common'
                                      '-issues,-and-recommendations#42-hits-where-the-spades-contig-contains'
                                      '-frameshifts\n. Default is: %(default)s.')
    parser_assemble.add_argument('--exonerate_skip_hits_with_internal_stop_codons',
                                 action='store_true',
                                 dest='skip_internal_stops',
                                 default=False,
                                 help='Skip Exonerate hits where the SPAdes sequence contains an internal in-frame '
                                      'stop codon. '
                                      'See:\nhttps://github.com/mossmatters/HybPiper/wiki/Troubleshooting,'
                                      '-common-issues,-and-recommendations#31-sequences-containing-stop-codons.\n'
                                      'A single terminal stop codon is allowed, but see option '
                                      '"--exonerate_skip_hits_with_terminal_stop_codons" below. Default is: '
                                      '%(default)s.')
    parser_assemble.add_argument('--exonerate_skip_hits_with_terminal_stop_codons',
                                 action='store_true',
                                 dest='skip_terminal_stops',
                                 default=False,
                                 help='Skip Exonerate hits where the SPAdes sequence contains a single terminal stop '
                                      'codon. Only applies when option '
                                      '"--exonerate_skip_hits_with_internal_stop_codons" is also provided. Only use '
                                      'this flag if your target file exclusively contains protein-coding genes with no '
                                      'stop codons included, and you would like to prevent any in-frame stop codons in '
                                      'the output sequences. Default is: %(default)s.')
    parser_assemble.add_argument('--merged',
                                 action='store_true',
                                 default=False,
                                 help='For assembly with both merged and unmerged (interleaved) reads. Default is: '
                                      '%(default)s.')
    parser_assemble.add_argument('--no_intronerate',
                                 action='store_true',
                                 dest='no_intronerate',
                                 default=False,
                                 help='Do not run intronerate to recover fasta files for supercontigs with introns (if '
                                      'present), and introns-only. Default is: %(default)s.')
    parser_assemble.add_argument('--keep_intermediate_files',
                                 action='store_true',
                                 dest='keep_intermediate_files',
                                 default=False,
                                 help='Keep all intermediate files and logs, which can be useful for '
                                      'debugging. Default action is to delete them, which greatly reduces the total '
                                      'file number. Default is: %(default)s.')
    parser_assemble.add_argument('--no_padding_supercontigs',
                                 action='store_true',
                                 dest='no_padding_supercontigs',
                                 default=False,
                                 help='If Intronerate is run, and a supercontig is created by concatenating multiple '
                                      'SPAdes contigs, do not add 10 "N" characters between contig joins. By default, '
                                      'Ns will be added. Default is: %(default)s.')
    parser_assemble.add_argument('--verbose_logging',
                                 action='store_true',
                                 dest='verbose_logging',
                                 default=False,
                                 help='If supplied, enable verbose logging. NOTE: this can increase the size of the '
                                      'log files by an order of magnitude. Default is: %(default)s.')
    parser_assemble.add_argument('--hybpiper_output', '-o',
                                 dest='output_folder',
                                 default=None,
                                 help='Folder for HybPiper output. Default is %(default)s.')
    parser_assemble.add_argument('--no_spades_eta',
                                 action='store_true',
                                 default=False,
                                 help='When SPAdes is run concurrently using GNU parallel, the "--eta" flag can '
                                      'result in many "sh: /dev/tty: Device not configured" errors written to stderr. '
                                      'Using this option removes the "--eta" flag to GNU parallel. Default is '
                                      '%(default)s.')
    parser_assemble.add_argument('--run_profiler',
                                 action='store_true',
                                 dest='run_profiler',
                                 default=False,
                                 help='If supplied, run the subcommand using cProfile. Saves a *.csv file of results. '
                                      'Default is: %(default)s.')

    # Set defaults for subparser <parser_assemble>:
    parser_assemble.set_defaults(blast=True)

    return parser_assemble


def add_stats_parser(subparsers):
    """
    Parser for hybpiper_stats, which now includes running get_seq_lengths

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_stats = subparsers.add_parser('stats', help='Gather statistics about the HybPiper run(s)')
    group_1 = parser_stats.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna',
                         dest='targetfile_dna',
                         default=False,
                         help='FASTA file containing DNA target sequences for each gene. The fasta headers must '
                              'follow the naming convention: >TaxonID-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa',
                         dest='targetfile_aa',
                         default=False,
                         help='FASTA file containing amino-acid target sequences for each gene. The fasta headers '
                              'must follow the naming convention: >TaxonID-geneName')
    parser_stats.add_argument('sequence_type',
                              choices=['gene', 'GENE', 'supercontig', 'SUPERCONTIG'],
                              help='Sequence type (gene or supercontig) to recover lengths for.')
    parser_stats.add_argument('namelist',
                              help="Text file with names of HybPiper output directories, one per line.")
    parser_stats.add_argument('--seq_lengths_filename',
                              default='seq_lengths',
                              help='File name for the sequence lengths *.tsv file. Default is <seq_lengths.tsv>.')
    parser_stats.add_argument('--stats_filename',
                              default='hybpiper_stats',
                              help='File name for the stats *.tsv file. Default is: <hybpiper_stats.tsv>')
    parser_stats.add_argument('--run_profiler',
                              action='store_true',
                              dest='run_profiler',
                              default=False,
                              help='If supplied, run the subcommand using cProfile. Saves a *.csv file of results')

    return parser_stats


def add_retrieve_sequences_parser(subparsers):
    """
    Parser for retrieve_sequences

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_retrieve_sequences = subparsers.add_parser('retrieve_sequences',
                                                      help='Retrieve sequences generated by the HybPiper run(s)')
    group_1 = parser_retrieve_sequences.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna',
                         dest='targetfile_dna',
                         default=False,
                         help='FASTA file containing DNA target sequences for each gene. The fasta headers must '
                              'follow the naming convention: >TaxonID-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa',
                         dest='targetfile_aa',
                         default=False,
                         help='FASTA file containing amino-acid target sequences for each gene. The fasta headers '
                              'must follow the naming convention: >TaxonID-geneName')
    parser_retrieve_sequences.add_argument('--sample_names',
                                           default=None,
                                           help='Directory containing Hybpiper output OR a file containing HybPiper '
                                                'output names, one per line.')
    parser_retrieve_sequences.add_argument('--single_sample_name',
                                           default=None,
                                           help='A single sample name to recover sequences for.')
    parser_retrieve_sequences.add_argument('sequence_type',
                                           choices=['dna', 'aa', 'intron', 'supercontig'],
                                           help='Type of sequence to extract.')
    parser_retrieve_sequences.add_argument('--hybpiper_dir',
                                           default=None,
                                           help='Specify directory containing HybPiper output.')
    parser_retrieve_sequences.add_argument('--fasta_dir',
                                           default=None,
                                           help='Specify directory for output FASTA files.')
    parser_retrieve_sequences.add_argument('--skip_chimeric_genes',
                                           action='store_true',
                                           dest='skip_chimeric',
                                           default=False,
                                           help='Do not recover sequences for putative chimeric genes. This only has '
                                                'an effect for a given sample if the option '
                                                '"--chimeric_stitched_contig_check" was provided to command '
                                                '"hybpiper assemble".')
    parser_retrieve_sequences.add_argument('--stats_file',
                                           default=None,
                                           help='Stats file produced by "hybpiper stats", required for selective '
                                                'filtering of retrieved sequences.')
    parser_retrieve_sequences.add_argument('--filter_by',
                                           action='append',
                                           nargs=3,
                                           metavar=('column', 'comparison', 'threshold'),
                                           help='Provide three space-separated arguments: 1) column of the stats_file '
                                                'to filter by, 2) "greater" or "smaller", 3) a threshold - either an '
                                                'integer (raw number of genes) or float (percentage of genes in '
                                                'analysis). This parameter can be supplied more than once to filter '
                                                'by multiple criteria.')
    parser_retrieve_sequences.add_argument('--run_profiler',
                                           action='store_true',
                                           dest='run_profiler',
                                           default=False,
                                           help='If supplied, run the subcommand using cProfile. Saves a *.csv file '
                                                'of results.')

    return parser_retrieve_sequences


def add_paralog_retriever_parser(subparsers):
    """
    Parser for paralog_retriever

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_paralog_retriever = subparsers.add_parser('paralog_retriever',
                                                     help='Retrieve paralog sequences generated by the HybPiper run(s)')
    parser_paralog_retriever.add_argument('namelist',
                                          help='Text file containing list of HybPiper output directories, '
                                               'one per line.')
    group_1 = parser_paralog_retriever.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna',
                         dest='targetfile_dna',
                         default=False,
                         help='FASTA file containing DNA target sequences for each gene. Used to extract unique gene '
                              'names for paralog recovery. The fasta headers must follow the naming convention: '
                              '>TaxonID-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa',
                         dest='targetfile_aa',
                         default=False,
                         help='FASTA file containing amino-acid target sequences for each gene. Used to extract '
                              'unique gene names for paralog recovery. The fasta headers must follow the naming '
                              'convention: >TaxonID-geneName')
    parser_paralog_retriever.add_argument('--fasta_dir_all',
                                          default='paralogs_all',
                                          help='Specify directory for output FASTA files (ALL). Default is: '
                                               '%(default)s.')
    parser_paralog_retriever.add_argument('--fasta_dir_no_chimeras',
                                          default='paralogs_no_chimeras',
                                          help='Specify directory for output FASTA files (no putative chimeric '
                                               'sequences). Default is: %(default)s.')
    parser_paralog_retriever.add_argument('--paralog_report_filename',
                                          default='paralog_report',
                                          help='Specify the filename for the paralog *.tsv report table. Default is: '
                                               '%(default)s.')
    parser_paralog_retriever.add_argument('--paralogs_above_threshold_report_filename',
                                          default='paralogs_above_threshold_report',
                                          help='Specify the filename for the *.txt list of genes with paralogs in '
                                               '<paralogs_list_threshold_percentage> number of samples. Default is: '
                                               '%(default)s.')
    parser_paralog_retriever.add_argument('--paralogs_list_threshold_percentage',
                                          type=float,
                                          default=0.0,
                                          help='Percent of total number of samples and genes that must have paralog '
                                               'warnings to be reported in the <genes_with_paralogs.txt> report file. '
                                               'The default is %(default)s, meaning that all genes and samples with '
                                               'at least one paralog warning will be reported')
    parser_paralog_retriever.add_argument('--no_heatmap',
                                          action='store_true',
                                          default=False,
                                          help='If supplied, do not create a heatmap figure. Default is: %(default)s.')
    parser_paralog_retriever.add_argument('--heatmap_filename',
                                          default='paralog_heatmap',
                                          help='Filename for the output heatmap, saved by default as a *.png file. '
                                               'Default is: %(default)s.')
    parser_paralog_retriever.add_argument('--figure_length',
                                          type=int,
                                          default=None,
                                          help='Length dimension (in inches) for the output heatmap file. Default is '
                                               'automatically calculated based on the number of genes.')
    parser_paralog_retriever.add_argument('--figure_height',
                                          type=int,
                                          default=None,
                                          help='Height dimension (in inches) for the output heatmap file. Default is '
                                               'automatically calculated based on the number of samples.')
    parser_paralog_retriever.add_argument('--sample_text_size',
                                          type=int,
                                          default=None,
                                          help='Size (in points) for the sample text labels in the output heatmap '
                                               'file. Default is automatically calculated based on the number of '
                                               'samples.')
    parser_paralog_retriever.add_argument('--gene_text_size',
                                          type=int,
                                          default=None,
                                          help='Size (in points) for the gene text labels in the output heatmap file. '
                                               'Default is automatically calculated based on the number of genes.')
    parser_paralog_retriever.add_argument('--heatmap_filetype',
                                          choices=['png', 'pdf', 'eps', 'tiff', 'svg'],
                                          default='png',
                                          help='File type to save the output heatmap image as. Default is: '
                                               '%(default)s.')
    parser_paralog_retriever.add_argument('--heatmap_dpi',
                                          type=int,
                                          default=100,
                                          help='Dots per inch (DPI) for the output heatmap image. Default is: '
                                               '%(default)s.')
    parser_paralog_retriever.add_argument('--no_xlabels',
                                          action='store_true',
                                          default=False,
                                          help='If supplied, do not render labels for x-axis (loci) in the saved '
                                               'heatmap figure. Default is: %(default)s.')
    parser_paralog_retriever.add_argument('--no_ylabels',
                                          action='store_true',
                                          default=False,
                                          help='If supplied, do not render labels for y-axis (samples) in the '
                                               'saved heatmap figure. Default is: %(default)s.')
    parser_paralog_retriever.add_argument('--run_profiler',
                                          action='store_true',
                                          dest='run_profiler',
                                          default=False,
                                          help='If supplied, run the subcommand using cProfile. Saves a *.csv file '
                                               'of results. Default is: %(default)s.')

    return parser_paralog_retriever


def add_gene_recovery_heatmap_parser(subparsers):
    """
    Parser for gene_recovery_heatmap

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_gene_recovery_heatmap = subparsers.add_parser('recovery_heatmap',
                                                         help='Create a gene recovery heatmap for the HybPiper run(s)')
    parser_gene_recovery_heatmap.add_argument('seq_lengths_file',
                                              help="Filename for the seq_lengths file (output of the 'hybpiper "
                                                   "stats' command)")
    parser_gene_recovery_heatmap.add_argument('--heatmap_filename',
                                              default='recovery_heatmap',
                                              help='Filename for the output heatmap, saved by default as a *.png file. '
                                                   'Defaults to "recovery_heatmap"')
    parser_gene_recovery_heatmap.add_argument('--figure_length',
                                              type=int,
                                              default=None,
                                              help='Length dimension (in inches) for the output heatmap file. '
                                                   'Default is automatically calculated based on the number of '
                                                   'genes')
    parser_gene_recovery_heatmap.add_argument('--figure_height',
                                              type=int,
                                              default=None,
                                              help='Height dimension (in inches) for the output heatmap file. '
                                                   'Default is automatically calculated based on the number of '
                                                   'samples')
    parser_gene_recovery_heatmap.add_argument('--sample_text_size',
                                              type=int,
                                              default=None,
                                              help='Size (in points) for the sample text labels in the output heatmap '
                                                   'file. Default is automatically calculated based on the '
                                                   'number of samples')
    parser_gene_recovery_heatmap.add_argument('--gene_text_size',
                                              type=int,
                                              default=None,
                                              help='Size (in points) for the gene text labels in the output heatmap '
                                                   'file. Default is automatically calculated based on the '
                                                   'number of genes')
    parser_gene_recovery_heatmap.add_argument('--heatmap_filetype',
                                              choices=['png', 'pdf', 'eps', 'tiff', 'svg'],
                                              default='png',
                                              help='File type to save the output heatmap image as. Default is *.png')
    parser_gene_recovery_heatmap.add_argument('--heatmap_dpi',
                                              type=int,
                                              default=100,
                                              help='Dot per inch (DPI) for the output heatmap image. Default is '
                                                   '%(default)d')
    parser_gene_recovery_heatmap.add_argument('--no_xlabels',
                                              action='store_true',
                                              default=False,
                                              help='If supplied, do not render labels for x-axis (loci) in the saved '
                                                   'heatmap figure')
    parser_gene_recovery_heatmap.add_argument('--no_ylabels',
                                              action='store_true',
                                              default=False,
                                              help='If supplied, do not render labels for y-axis (samples) in the '
                                                   'saved heatmap figure')
    parser_gene_recovery_heatmap.add_argument('--run_profiler',
                                              action='store_true',
                                              dest='run_profiler',
                                              default=False,
                                              help='If supplied, run the subcommand using cProfile. Saves a *.csv file '
                                                   'of results')

    return parser_gene_recovery_heatmap


def add_check_dependencies_parser(subparsers):
    """
    Parser for check_dependencies

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_check_dependencies = subparsers.add_parser('check_dependencies',
                                                      help='Run a check for all pipeline dependencies and exit')

    parser_check_dependencies.add_argument('--run_profiler',
                                           help='If supplied, run the subcommand using cProfile. Saves a *.csv file '
                                                'of results',
                                           action='store_true',
                                           dest='run_profiler',
                                           default=False)

    # Set defaults for subparser <check_dependencies>:
    parser_check_dependencies.set_defaults(logger=None)

    return parser_check_dependencies


def add_check_targetfile_parser(subparsers):
    """
    Parser for check_targetfile

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_check_target_file = subparsers.add_parser('check_targetfile',
                                                     help='Check the target file for issues, then exit')
    group_1 = parser_check_target_file.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna',
                         dest='targetfile_dna',
                         default=False,
                         help='FASTA file containing DNA target sequences for each gene. The fasta headers must '
                              'follow the naming convention: >TaxonID-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa',
                         dest='targetfile_aa',
                         default=False,
                         help='FASTA file containing amino-acid target sequences for each gene. The fasta headers '
                              'must follow the naming convention: >TaxonID-geneName')
    parser_check_target_file.add_argument('--no_terminal_stop_codons',
                                          action='store_true',
                                          default=False,
                                          help='When testing for open reading frames, do not allow a translated frame '
                                               'to have a single stop codon at the C-terminus of the translated '
                                               'protein sequence. Default is False.')
    parser_check_target_file.add_argument('--sliding_window_size',
                                          type=int,
                                          default=None,
                                          help='Number of characters (single-letter DNA or amino-acid codes) to '
                                               'include in the sliding window when checking for sequences with '
                                               'low-complexity-regions.')
    parser_check_target_file.add_argument('--complexity_minimum_threshold',
                                          type=float,
                                          default=None,
                                          help='Minimum threshold value. Beneath this value, the sequence in the '
                                               'sliding window is flagged as low complexity, and the corresponding '
                                               'target file sequence is reported as having low-complexity regions.')
    parser_check_target_file.add_argument('--run_profiler',
                                          help='If supplied, run the subcommand using cProfile. Saves a *.csv file '
                                               'of results',
                                          action='store_true',
                                          dest='run_profiler',
                                          default=False)

    # Set defaults for subparser <check_target_file>:
    parser_check_target_file.set_defaults(logger=None)

    return parser_check_target_file


def add_fix_targetfile_parser(subparsers):
    """
    Parser for fix_targetfile

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_fix_target_file = subparsers.add_parser('fix_targetfile',
                                                   help='Fix and filter the target file, then exit')
    group_1 = parser_fix_target_file.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna',
                         dest='targetfile_dna',
                         default=False,
                         help='FASTA file containing DNA target sequences for each gene. The fasta headers must '
                              'follow the naming convention: >TaxonID-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa',
                         dest='targetfile_aa',
                         default=False,
                         help='FASTA file containing amino-acid target sequences for each gene. The fasta headers '
                              'must follow the naming convention: >TaxonID-geneName')
    parser_fix_target_file.add_argument('control_file',
                                        help='The *.ctl file, as output by the command "hybpiper check_targetfile".')
    parser_fix_target_file.add_argument('--no_terminal_stop_codons',
                                        action='store_true',
                                        default=False,
                                        help='When testing for open reading frames, do not allow a translated frame '
                                             'to have a single stop codon at the C-terminus of the translated '
                                             'protein sequence. Default is False. If supplied, this parameter will '
                                             'override the setting in the *.ctl file.')
    parser_fix_target_file.add_argument('--allow_gene_removal',
                                        action='store_true',
                                        default=False,
                                        help='Allow frame-correction and filtering steps to remove all representative '
                                             'sequences for a given gene. Default is False; HybPiper will exit with an '
                                             'information message instead. If supplied, this parameter will '
                                             'override the setting in the *.ctl file.')
    parser_fix_target_file.add_argument('--reference_protein_file',
                                        default=None,
                                        help='If a given DNA sequence can be translated in more than one forward frame '
                                             'without stop codons, choose the translation that best matches the '
                                             'corresponding reference protein provided in this fasta file. The fasta '
                                             'headers must follow the naming convention: >TaxonID-geneName')
    parser_fix_target_file.add_argument('--maximum_distance',
                                        default=0.5,
                                        type=utils.restricted_float,
                                        metavar='FLOAT',
                                        help='When comparing candidate DNA translation frames to a reference protein, '
                                             'the maximum distance allowed between the translated frame and the '
                                             'reference sequence for any candidate translation frame to be selected. '
                                             'Useful to filter out sequences with frameshifts that do NOT introduce '
                                             'stop codons. 0.0 means identical sequences, 1.0 means completely '
                                             'different sequences. Default is 0.5')
    parser_fix_target_file.add_argument('--filter_by_length_percentage',
                                        default=0.0,
                                        type=utils.restricted_float,
                                        metavar='FLOAT',
                                        help='If more than one representative sequence is present for a given gene, '
                                             'filter out sequences shorter than this percentage of the longest gene '
                                             'sequence length. Default is 0.0 (all sequences retained).')
    parser_fix_target_file.add_argument('--keep_low_complexity_sequences',
                                        action='store_true',
                                        default=False,
                                        help='Keep sequences that contain regions of low complexity, as identified by '
                                             'the command "hybpiper check_targetfile". Default is to remove these '
                                             'sequences.')
    parser_fix_target_file.add_argument('--alignments',
                                        action='store_true',
                                        default=False,
                                        help='Create per-gene alignments from the final fixed/filtered target file '
                                             'sequences. Note that DNA sequences will be translated prior to '
                                             'alignment.')
    parser_fix_target_file.add_argument('--concurrent_alignments',
                                        default=1,
                                        type=int,
                                        metavar='INTEGER',
                                        help='Number of alignments to run concurrently. Default is 1.')
    parser_fix_target_file.add_argument('--threads_per_concurrent_alignment',
                                        default=1,
                                        type=int,
                                        metavar='INTEGER',
                                        help='Number of threads to run each concurrent alignment with. Default is 1.')
    parser_fix_target_file.add_argument('--write_all_fasta_files',
                                        default=False,
                                        action='store_true',
                                        help='If provided, *.fasta files will be written for sequences removed from '
                                             'the fixed/filtered target file, according to filtering categories '
                                             '(length threshold, low-complexity regions, etc.). By default, '
                                             'these files will not be written.')
    parser_fix_target_file.add_argument('--verbose_logging',
                                        help='If supplied, enable verbose logging. NOTE: this will increase the size '
                                             'of the log files.',
                                        action='store_true',
                                        dest='verbose_logging',
                                        default=False)
    parser_fix_target_file.add_argument('--run_profiler',
                                        help='If supplied, run the subcommand using cProfile. Saves a *.csv file '
                                             'of results',
                                        action='store_true',
                                        dest='run_profiler',
                                        default=False)

    # Set defaults for subparser <fix_target_file>:
    parser_fix_target_file.set_defaults(logger=None)

    return parser_fix_target_file


def add_filter_by_length_parser(subparsers):
    """
    Parser for filter_by_length

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_filter_by_length = subparsers.add_parser('filter_by_length',
                                                    help='Filter the sequences output by command "hybpiper '
                                                         'retrieve_sequences" by length (absolute and relative to mean)')

    parser_filter_by_length.add_argument('sequence_type',
                                         choices=["dna", "aa", "supercontig", "intron"],
                                         help='File sequence type for all FASTA files to filter in current directory. '
                                              'For example, the amino-acid output of HybPiper would be: aa')
    group_1 = parser_filter_by_length.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--denylist',
                         default=None,
                         help='Text file containing gene-sample combinations to omit.\nThe format of the file should '
                              'be one gene per line, a tab, and then a comma-delimited list of samples to disallow: '
                              '\n\n\tgene[tab]sample,sample,sample ')
    group_1.add_argument('--seq_lengths_file',
                         help='Filename for the seq_lengths file (output of the "hybpiper stats" command), with a list '
                              'of genes in the first row, mean target lengths in the second row, and sample recovery '
                              'in other rows.')
    parser_filter_by_length.add_argument('--denylist_filename',
                                         default='denylist.txt',
                                         type=str,
                                         help='File name for the "deny list" text file (if written). Default is '
                                              '<denylist.txt>')
    parser_filter_by_length.add_argument('--length_filter',
                                         default=0,
                                         type=int,
                                         help='Minimum length to allow a sequence in nucleotides for DNA or amino '
                                              'acids for protein sequences')
    parser_filter_by_length.add_argument('--percent_filter',
                                         default=0,
                                         type=float,
                                         help='Minimum fraction (between 0 and 1) of the mean target length to allow '
                                              'a sequence for a gene. Lengths taken from HybPiper stats file.')
    parser_filter_by_length.add_argument('--sequence_dir',
                                         default=None,
                                         help='Specify directory containing sequences output by the "hybpiper '
                                              'retrieve_sequences" command. Default is to search in the current '
                                              'working directory')
    parser_filter_by_length.add_argument('--filtered_dir',
                                         default=None,
                                         help='Specify directory for output filtered FASTA files. Default is to write '
                                              'to the current working directory')
    parser_filter_by_length.add_argument('--run_profiler',
                                         help='If supplied, run the subcommand using cProfile. Saves a *.csv file '
                                              'of results',
                                         action='store_true',
                                         dest='run_profiler',
                                         default=False)

    # Set defaults for subparser <filter_by_length>:
    parser_filter_by_length.set_defaults(logger=None)

    return parser_filter_by_length
