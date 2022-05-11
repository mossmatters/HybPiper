# Changelog

**2.0 release candidate** *March, 2022*

This update involves a substantial refactor of the HybPiper pipeline, with changes to the internal code, additional functionality, and additional output. Changes include:

***New dependencies***: 

- Python 3.6 or later
- seaborn (Python library)
- matplotlib (Python library)
- progressbar2 (Python library)
- BioPython 1.80 (Python library).
- pandas (Python library)
- DIAMOND
- BBtools (BBmap.sh, BBmerge.sh)

***Code refactor***:

- HybPiper now has an installation process (conda or pip via setup.py), rather than scripts being called directly.
- After installation, all HybPiper commands are accessed via the main command `hybpiper`, followed by a subcommand e.g. `assemble`.
- The`reads_first.py` module has been renamed to `assemble.py`, and now imports other pipeline modules rather than calling them as external Python scripts.
- The `exonerate_hits.py` module has been substantially re-written to use the BioPython SearchIO Exonerate text parser, as this allows much more data recovered from Exonerate search results. 
- The `intronerate.py` module has been removed; this functionality has been moved to `exonerate_hits.py`.
- The `paralog_investigator.py` module has been removed; this functionality has been moved to `exonerate_hits.py`.
- The `gene_recovery_heatmap_ggplot.R` and `gene_recovery_heatmap.R` R scripts have been replaced by the Python module `gene_recovery_heatmap.py`.
- The `get_seq_lengths.py` script has been removed; this functionality is now performed by the module `hybpiper_stats.py` (command `hybpiper stats`). The file `seq_lengths.tsv` is written by default; this filename can be changed via the parameter `--seq_lengths_filename`.  

***Additional features/functionality***:

- The program DIAMOND can be used in place on BLASTX when mapping reads to target sequences.
- When using BLAST or DIAMOND, the `hybpiper stats` command now calculates the enrichment efficiency; previously this was only calculated when using BWA.
- The `hybpiper assemble` commands can now accept read files in compressed gzip format (suffix `*.gz`).
- Logging when running `hybpiper assemble` has been unified and extended to provide additional debugging information. A single log file is written per-sample in the sample directory e.g. `EG30_reads_first_2021-12-02-10_45_56.log`. 
- Checks for all dependencies are now run by default when `hybpiper assemble` is run.
- The `hybpiper assemble` command checks that the provided target file is formatted correctly and can be translated as expected (in the case of a nucleotide target file). Any issues are printed to screen and details are logged to file.
- All Exonerate searches are now performed with the option `--refine full`; in the case of failure, a fallback run without this parameter is performed.
- In some situations HybPiper creates a gene sequence by 'stitching' together Exonerate hits from different SPAdes contigs. In some scenarios this can result in hit from different paralogs being joined together. HybPiper now performs a test (work in progress, currently sensitivity is low) to search for such 'chimeric' sequences, and provides warnings in the file `{sample_name}_genes_derived_from_putative_chimera_stitched_contig.csv` within each sample folder. Note that this chimera test is only performed in cases where a stitched contig has been created from multiple contigs, and paired-end reads are provided.
- By default, the SPAdes assembly folder is now deleted for each gene after contigs have been recovered. The user no longer needs to run `cleanup.py` after each run, and this script has been removed. Deleting the SPAdes directory dramatically reduces the total number of files produced by a completed run of HybPiper, which can be very useful when running it on and HPC with file number limits. To retain the SPAdes directory (i.e. for debugging purposes), the flag `--keep_spades_folder` can be used.
- When running Intronerate (via flag `-run_intronerate` with command `hybpiper assemble`, a block of 10 'N' characters is inserted into supercontigs at locations where different SPAdes contigs have been concatenated. This behaviour can be turned off via the flag `--no_padding_supercontigs`.
- In cases where HybPiper recovers sequence for multiple non-contiguous segments of a gene, the gaps between the  segments will be padded with a number of 'N' characters. The number of Ns corresponds to the number of amino acids in 'best' protein reference for that gene that do not have corresponding SPAdes contig hits, multiplied by 3 to convert to nucleotides.
- The command `hybpiper stats` now writes pipeline run statistics directly to file, rather than to standard out.
- The command `hybpiper paralog_retriever` now writes a table of gene-vs-sample with a matrix of paralog counts, and produces a heatmap image file from this table. A text report file is also produced, listing gene and sample names that contain paralogs in above a given threshold percentage. 
  
  
***New options/flags***:

For subcommand **`hybpiper assemble`**:
- `--diamond`. When provided, DIAMOND (with default sensitivity `fast`) will be used to map reads against the bait/target file, rather than the default BLASTX.
- `--diamond_sensitivity`. DIAMOND will be run with the provided sensitivity (options are `mid-sensitive`, `sensitive`, `more-sensitive`, `very-sensitive`, `ultra-sensitive`).
- `--run_intronerate`. When used, the function `intronerate()` in module `exonerate_hits.py` will be run to recover introns and supercontigs where possible.
- `--merged`. When used, R1 and R2 reads will be merged using `BBmerge.sh`, where possible. Both the merged and the remaining unmerged reads will be used for SPAdes assembly.
- `--no_stitched_contigs`. When used, gene sequences will comprise the longest Exonerate hit from a single SPAdes contig; no contig/hit stitching will be attempted.
- `--paralog_min_length_percentage`. Corresponds to the minimum percentage length for a SPADes contig Exonerate hit (vs the reference query sequence length) for it to be flagged and recovered as a 'long' paralog. Previously this parameter was hardcoded to 0.75. This parameter effected both types of paralog warnings (by length, and by depth).
- `--bbmap_subfilter`. Ban BBmap alignments with more than this many substitutions when searching for chimeric stitched contigs. Default is 7.
- `--bbmap_threads`. The number of threads to use for BBmap when searching for chimeric stitched contigs. Default is 2.
- `--bbmap_memory`. The amount of memory (RAM) in GB to use for BBmap when searching for chimeric stitched contigs. Default is 1.
- `--chimeric_stitched_contig_edit_distance`. Minimum number of differences between one read of a read pair vs a stitched contig reference for a read pair to be flagged as discordant. Default is 5. 
- `--chimeric_stitched_contig_discordant_reads_cutoff`. Minimum number of discordant reads pairs required to flag a stitched contig as a potential chimera of contigs from multiple paralogs. Default is 5.
- `--keep_spades_folder`. If used, the SPAdes assembly folder for each gene will not be deleted. Note that previous versions of HybPiper retained the SPAdes assembly folders by default; they could previously be removed after running `reads_first.py` using the `cleanup.py` script.  
- `--keep_exonerate_logs.`. If used, retain the `exonerate_hits.py` module log within each gene folder. Default behaviour is to delete each log after it has been copied to the main log file in the sample directory.

For subcommand **`hybpiper retrieve_sequences`**:

- `--single_sample_name`. Allows the user to specify a single sample name; sequences from this sample only will  be recovered in single `*.fasta` file.
- `--skip_chimeric_genes`. If used, skip retrieval of any genes/sample for which the stitiched contig was flagged as chimeric (potentially derived from multiple paralogs).

For subcommand **`hybpiper stats`**:

- A positional parameter for the `sequence type` has been added, with the options {gene, supercontig}. Note that all length statistics are now provived as number of nucleotides, regardless of whether an amino-acid or nucleotide bait file was used.
- `--stats_filename`. The statistics are now written directly to file. This parameter can be used to specify the filename; default is `hybpiper_stats.tsv`.


For subcommand **`hybpiper recovery_heatmap`**:

- `seq_lengths_file`. A positional parameter providing the filename for the `seq_length.tsv` (default) output of `hybpiper stats` 
- `--heatmap_filename`. Supply a filename for the heatmap image file. Default is `heatmap`.
- `--figure_length`. Supply a length dimension (in inches) for the heatmap image file. Default is auto-calculated.
- `--figure_height`. Supply a height dimension (in inches) for the heatmap image file. Default is auto-calculated.
- `--sample_text_size`. Supply a text size (in points) for the heatmap image file sample names. Default is auto-calculated.
- `--gene_text_size`. Supply a text size (in points) for the heatmap image file gene names. Default is auto-calculated.
- `--heatmap_filetype`. Provide a filetype for the heatmap image file {png,pdf,eps,tiff,svg}. Default is `png`.
- `--heatmap_dpi`. Supply a Dots Per Inch value (in DPI) for the heatmap image file. Default is 300.


For subcommand **`hybpiper paralog_retriever`**:

- `--fasta_dir_all`. Supply a directory name for the paralog `*.fasta` files when recovering all sequences included putative chimeric stitched contig sequences (the latter are recovered when paralogs are not present). Default is `paralogs_all`.
- `--fasta_dir_no_chimeras`. Supply a directory name for the paralog `*.fasta` files when recovering all sequences _except_ for putative chimeric stitched contig sequences. Default is `paralogs_no_chimeras`.
- `--paralogs_above_threshold_report_filename`. Specify the filename for a text report file listing genes and samples with paralogs in greater than a given percentage threshold of genes/samples. Default is `genes_with_paralogs.txt`.
- `--paralogs_list_threshold_percentage`. Percent of total number of samples and genes that must have paralog warnings to be reported in the `genes_with_paralogs.txt` (default) file.
- `--heatmap_filename`. Supply a filename for the heatmap image file. Default is `paralog_heatmap`.
- `--figure_length`. Supply a length dimension (in inches) for the heatmap image file. Default is auto-calculated.
- `--figure_height`. Supply a height dimension (in inches) for the heatmap image file. Default is auto-calculated.
- `--sample_text_size`. Supply a text size (in points) for the heatmap image file sample names. Default is auto-calculated.
- `--gene_text_size`. Supply a text size (in points) for the heatmap image file gene names. Default is auto-calculated.
- `--heatmap_filetype`. Provide a filetype for the heatmap image file {png,pdf,eps,tiff,svg}. Default is `png`.
- `--heatmap_dpi`. Supply a Dots Per Inch value (in DPI) for the heatmap image file. Default is 300.

***The following options/flags have been changed or removed***:

- The parameter `--baitfile` or `-b` has been replaced with `--targetfile` or `-t`, to more accurately describe the contents of the file provided. 
- The parameter `--length_pct` has been removed as a parameter to `hybpiper assemble` and is no longer used in internal code.
- The parameter `--thresh` (percent identity threshold for retaining Exonerate hits) for `hybpiper assemble` now defaults to 55 (previously 65).
- The flag `--check-depend` has been replaced by `--check_dependencies_only` for command `hybpiper assemble`. Dependency checking is now also performed every time the `hybpiper assemble` command is run.

***The following output files have been changed or removed***:

- XXX

***The following output files/folders have been added***:

- The script `paralog_retreiver.py` now write paralogs to two folders - one with all sequences, and the other without putative chimeric sequences.
- The file containing paralog warnings (produced when multiple long contigs are present for a gene) has been renamed from `genes_with_paralog_warnings.txt` to `<sample_name>_genes_with_long_paralog_warnings.txt`. 
- In addition to the standard paralog warning produced when multiple long contigs are present for a gene, the `hybpiper assemble` command now provides a paralog warning when multiple short contigs are present which together cover the reference sequence for a given gene at a depth >1, across a given percentage length (default 75%) of the reference. These warning are written to each sample directory to the file `<sample_name>_genes_with_paralog_warnings_by_contig_depth.csv`
- `chimera_test_diagnostic_reads.sam`.
- XXX Other intronerate files, if we end up keeping them?
- XXX Heatmap and extra reports for paralog_retriever
- XXX retrieve sequences output folders

    
**1.3.2** *February, 2020*

- Fix for [Issue 41](https://github.com/mossmatters/HybPiper/issues/41) a problem in `intronerate.py` when attempting to resolve overlapping gene annotations.
- Add support in `retrieve_sequences.py` for a list of HybPiper outputs rather than everything in a directory 
- Add support for supercontigs in `get_seq_lengths.py`
- Remove integer requirement for `--cov_cutoff` to accommodate `auto` and `off` settings in Spades.

**1.3.1 IMPORTANT BUG FIX** *August, 2018*

- Reverts a change in 1.3 that changed the behavior when there is only one hit in `exonerate`. Strandedness was not handled properly, leading some sequences to be returned in reverse complement. **Sequences recovered using 1.3 should be re-run with `reads_first.py --no-blast --no-distribute --no-assemble` to repeat just the exonerate step.**


**1.3 The Herbarium Update** *January, 2018*

### Features 
- Added `--exclude`  flag to be the inverse of `--target`: all sequences with the specified string will not be used as targets for exon extraction (they will still be used for read-mapping). Useful if you want to add supercontig sequence to the target file, but not use it for exon extraction.

- Added `--addN` to `intronerate.py`. This feature will add 10 N characters in between joined contig when recovering the supercontig. This is useful for identifying where the intron recovery fails, and for annotation processing (i.e. for GenBank).

- Added a new version of the heatmap script, `gene_recovery_heatmap_ggplot.R`. This script is much simpler and produces nice color PNG images, but struggles a bit on PDF output. The original heatmap script is stil included. *Thanks to Paul Wolf for the ggplot code!*


### Bug Fixes
- Fixed misassembly of supercontigs when there are multiple alignments to different parts of the same exon.
- Fixed poor filtering of GFF results to produce intron/exon annotation.
- Fixed non-propogation of exonerate parameters



**1.2.1** *September, 2017*

### Bug Fixes
- Fixed assembly issue when gene does not have unpaired reads.
- Fixed distribution of targets when unpaired reads present.
- Fixed use of unpaired reads to detect best target.


**1.2** *May, 2017*

### Features 

- Added `--unpaired` flag. When using paired-end sequencing reads, a third read file may be specified with this flag. Reads will be mapped to targets separately, but will be used along with paired reads in contig assembly.

- Added `--target` flag. Adds the ability to choose which of the reference sequences is used for each gene. If `--target` is a file (tab-delimited file with one gene and one target name per line), HybPiper will use that. Otherwise `--target` can be the name of one reference. HybPiper will only use targets with the specified name in the Alignment/Exon Extraction phase. All other targets for that locus will only be used in the Mapping/Read Sorting phase.

- Added `--timeout` flag, which uses GNU Parallel to kill processes (i.e. Spades or Exonerate) if they take X percent longer than average. Use if there are a lot of stuck jobs (`--timeout 1000`)
- Python 3 compatibility

### Bug Fixes
- Can accommodate Solexa FASTQ paired headers
- Fixed `spades_runner.py` not recognizing `--cpu` on redos
- Prints more meaningful messages for some common errors
- Can accommodate `prefix` not being in current directory
- Deletes sorted reads on restart to prevent double counting reads.
- `spades_runner.py` will now respect `--kvals`
- Added initial call to log for `reads_first.py`

---
**1.1** *May, 2016*: Release associated with manuscript in *Applications in Plant Sciences*.

- Added `paralog_investigator.py`, which detects and extracts long exons from putative paralogs in all genes in one sample.

- Added `paralog_retriever.py`, which retrieves sequences generated by `paralog_investigator.py` for many samples (or the coding sequence generated by `exonerate_hits.py` if no paralogs are detected).

- Added a test_dataset of 13 genes for 9 samples, and a shell script to run the test data through the main script and several post-processing scripts.

- Fixed bug involving calling HybPiper with a relative path such as: `../reads_first.py`

- `reads_first.py --check_depend` now checks for SPAdes, BWA, and Samtools

- Full revision of README, which is now shorter. Full tutorials on installing and running HybPiper are now on the Wiki.



**1.0** *Feb, 2016*: Initial fully-featured release associated with submission of manuscript to *Applications in Plant Sciences*.

-- Sequence assembly now uses SPAdes rather than Velvet + CAP3