# Changelog

**2.1.8** *25th June, 2024*

- Add new subcommand `hybpiper filter_by_length`, used to filter the sequence output of `hybpiper retrieve sequences` by absolute length and/or length relative to mean length in target file representatives. This is done on a per-sample/per-gene basis, rather than the sample-level filtering available in `hybpiper retrieve_sequences`. See [wiki](https://github.com/mossmatters/HybPiper/wiki#hybpiper-filter_by_length) for more information.
- Update the regex used to check target file fasta header formatting, to capture scenarios where a name contains multiple dashes and also ends with a dash.
- In the `fix_targetfile.py` module, remove the import of `Bio.Align.Applications.MafftCommandline` and call `MAFFT` via `subprocess` (see [issue#147](https://github.com/mossmatters/HybPiper/issues/147)).
- In the `gene_recovery_heatmap.py` module, cast the dataframe from the `seq_lengths_file` to object `dtype` to avoid a [deprecation warning](https://pandas.pydata.org/docs/whatsnew/v2.1.0.html#deprecated-silent-upcasting-in-setitem-like-series-operations) .
- Add option `--no_heatmap` to command `hybpiper paralog_retriever` (see [issue#150](https://github.com/mossmatters/HybPiper/discussions/150)).
- Fix an Exonerate-related debug message in `exonerate_hits.py`.

**2.1.7** *09th May, 2024*

- The flag `--run_intronerate` was removed from the `hybpiper assemble` command in the `run_hybpiper_test_dataset.sh` file.
- Removed the legacy check and attempted download of the test dataset in the `run_hybpiper_test_dataset.sh` file.
- Added a check to `hybpiper stats` and `hybpiper retrieve_sequences` to ensure sample names in the `namelist.txt` file do not contain forward slashes [issue#143](https://github.com/mossmatters/HybPiper/issues/143).
- When checking for putative chimeric gene sequences in `hybpiper retrieve_sequences` and `hybpiper paralog_retriever`, generate a warning rather than an error if the file `<sample_name>_genes_derived_from_putative_chimeric_stitched_contig.csv` can't be found for a given sample. This file will not be written if no gene sequences were produced for this sample (i.e. no reads mapped, no SPAdes contigs, no sequences extracted from SPAdes contigs via Exonerate).
- Check that target file FASTA headers do not contain quotation marks (`"` or `'`); [issue#125](https://github.com/mossmatters/HybPiper/issues/125).
- Updated the installation instructions in the README and Wiki to use the Bioconda package, and added installation instruction for Macs with Apple Silicon (M1/M2/M3 chips).
- Fixed a bug in `exonerate_hits.py` that meant that hits were not always trimmed to start with the first amino-acid with full alignment identity. This bug could potentially have had an effect on output sequences only if the values for `--exonerate_hit_sliding_window_size` and/or `--exonerate_hit_sliding_window_thresh` were changed from default values.
- Use `importlib.metadata` rather than `pkg_resources` for module version checks, due to deprecation of the latter. 

**2.1.6* * *19th July, 2023*

- **Intronerate is now run by default**. The flag `--run_intronerate` for subcommand `hybpiper assemble` has been changed to `--no_intronerate`. 
- If Intronerate fails, failed genes and errors will be printed and logged; the exonerate_contigs step of the pipeline will continue.
- Updated error handling and logging for the exonerate_contigs step of the pipeline.
- Change default DPI of heatmaps to 100 (previously 150) for `hybpiper recovery_heatmap` and `hybpiper paralog_retriever`
- Enforce rendering of all loci (x-axis) and sample (y-axis) labels in heatmaps; previously, matplotlib/seaborn would dynamically drop labels if they were too closely spaced.
- Added flags `--no_xlabels` and `--no_ylabels` for `hybpiper recovery_heatmap` and `hybpiper paralog_retriever`; turns off rendering of the corresponding labels in the saved figures.
- If the auto-calculated size of heatmaps for `hybpiper recovery_heatmap` and `hybpiper paralog_retriever` is greater than the maximum number of pixels (65536) in either/or length and height, resize the figure to 400 inches and 100 DPI. Note that large datasets can fail to render fully in the saved figure even if the pixel dimensions are less than the maximum (see e.g. https://stackoverflow.com/questions/64393779/how-to-render-a-heatmap-for-a-large-array), but reducing the size/DPI further allows the full figure to be rendered. 
- Added module `version.py` for a single location of HybPiper version number.
- Print and log HybPiper version when calling all subcommands.
- Added column 'TotalBasesRecovered' to the `hybpiper stats` report, listing the total number of nucleotides recovered for each sample (not counting N characters). Added 'TotalBasesRecovered' as a filtering option in `hybpiper retrieve_sequences`.

**2.1.5** *21st June, 2023*

- Bugfix: fixed an issue in `exonerate_hits.py` that could result in initial Exonerate hits being trimmed too aggressively at their 3' ends.
- Bugfix: fixed an issue in `exonerate_hits.py` that could introduce minor insertions in to the supercontig (concatenated exon and partial intron) sequence used when running Intronerate.  

**2.1.4** *5th June, 2023*

- Bugfix: fixed an issue when using `--run_intronerate` that could cause an error and result in no `*.FNA` sequence being produced for some genes.  


**2.1.3** *23rd March, 2023*

- Log platform and ulimit details for debugging purposes.
- Added parameter `--exonerate_hit_sliding_window_thresh` to `hybpiper assemble`. This value (default is 55) is used as the similarity threshold within the sliding window when trimming the 5' and 3' termini of Exonerate hits (i.e., the filter added in version 2.1.2). Previously, this filter used the value from the `--thresh` parameter, which is also used to perform filtering of entire Exonerate hits based on global alignment similarity.
- Bugfix: in `exonerate_hits.py`, ensure that similarity-filtered Exonerate hits are ordered by query start THEN query end; in cases where hits have identical query start coordinates, sorting only by query start can return the hits in a different order from one run to another, sometimes causing issues when trimming overlaps and concatenating to create a stitched-contig.
- Add parameter `--hybpiper_output` / `-o` to `hybpiper assemble`. If a directory is supplied using this parameter, `hybpiper assemble` will create it and use it for all output.


**2.1.2** *31st January, 2023*

- Removed the call to `time` when running SPAdes assemblies via `parallel` (see https://github.com/mossmatters/HybPiper/issues/109).
- Added a filter/trimming step to `exonerate_hits.py`. This process applies a sliding window to Exonerate hit alignments that have already passed a global similarity threshold (parameter `--thresh` for command `hybpiper assemble`). The 5' and 3' termini of SPAdes contig hits are trimmed if the similarity within the sliding window is below the `--thresh` value. This removes putative spurious 5' and 3' hit sequence produced by Exonerate alignments, and results in more accurate output `*.FNA` and `*.FAA` sequences. The sliding window size can be adjusted using the `hybpiper assemble` parameter `--exonerate_hit_sliding_window_size`; default value is 3 (i.e., 3 amino acids / 9 nucleotides).
- HybPiper now provides a warning if any output sequence contains internal stop codons, and writes the corresponding gene names to the file `<prefix>_genes_with_non_terminal_stop_codons.txt`.
- Fixed a bug that meant reads with older format headers (i.e., suffixes `/1` or `/2`) weren't processed properly when using BLASTx/DIAMOND (see https://github.com/mossmatters/HybPiper/issues/108).


**2.1.1** *12th December, 2022*

- When mapping reads with DIAMOND via `hybpiper assemble --diamond`, remove the `gunzip` step and on-the-fly fastq to fasta conversion (as DIAMOND supports both `*.fastq` and `*.gz` input). Further, pass the value of the `hybpiper assemble` parameter `--cpu`  directly to the `--threads` parameter of the `diamond blastx` command; do not run `diamond` via GNU parallel. See issue #104. 


**2.1.0** *1st December, 2022*

- The subcommand `hybpiper check_targetfile` now writes a `*.ctl file`; see wiki for details.
- Added the new subcommand `hybpiper fix_targetfile`; see wiki for details.
- New dependency: MAFFT. Used for optional alignments when running `hybpiper fix_targetfile`; see wiki for details.
- The flag `---distribute_high_mem` has been changed to `--distribute_low_mem` for subcommand `hybpiper assemble`.  Default is off, i.e., reads are now distributed using the faster, more memory intensive approach.
- The `time` command has been removed from all BWA and BLAST command strings (see https://github.com/mossmatters/HybPiper/issues/89).
- Update minimum Python version from 3.6 to 3.7, to match minimum requirements for Biopython 1.80 release.


**2.0.3** *16th November, 2022*

- The calculation for automatically determining the default heatmap dimensions has been changed, to prevent label trimming when using large numbers of samples and genes (pull request from LPDagallier).
- The default DPI of the heatmap `*.png` file has been reduced from 300 to 150.


**2.0.2** *14th November, 2022*

- Added the `--run_profiler` flag to the `check_dependencies` subcommand (bugfix). 
- Moved to correct semantic versioning; removed the 'release candidate build x' string from version description, and bumped to version 2.0.2.


**2.0.1 release candidate build 13** *12th September, 2022*

- Added a `--run_profiler` option to all subcommands; when used, the given subcommand is run with cProfile and the data is saved to a *.csv file.
- If a DNA target file is provided but the flag `--bwa` is omitted, a translated target file is now written directly to the sample directory (rather than the same directory as the input DNA file) with the generic name `translated_target_file.fasta`. This prevents file overwrite/access issues when processing multiple samples concurrently (e.g. using Nextflow).
- Added a check for sequence fasta files that exist but are empty to `stats`, `paralog_retriever` and `retrieve_sequences`; a warning is printed if such a file is found.


**2.0.1 release candidate build 12** *6th July, 2022*

- When running `hybpiper stats`, empty lines in the namelist.txt file are ignored.


**2.0.1 release candidate build 11** *29th June, 2022*

- Adding hybpiper_dir search for single sample recovery when running `hybpiper retrieve_sequences`. Adding sample name to recovered fasta file.
-  Adding `--memory 1024` to `spades.py` command (macOS Monterey compatibility pre SPAdes version 3.15.4)


**2.0.1 release candidate build 8** *24th June, 2022*

- Capture errors relating to malformed input `*.fastq` files during read distribution, and log error to file.
- When running `hybpiper stats`, if no input read count file can be found, log issue and sample to screen along with instructions, and exit.


**2.0.1 release candidate build 8** *16th June, 2022*

- If mapping is performed with BWA and the BWA mapping step fails, remove any `*.bam` file that has been produced.


**2.0 release candidate** *March, 2022*

This update involves a substantial refactor of the HybPiper pipeline, with changes to the internal code, additional functionality, and additional output. Changes include:


# New and updated dependencies
* [Python](https://www.python.org/downloads/) 3.6 or later, along with the Python libraries:
    * [seaborn](https://seaborn.pydata.org/installing.html)
    * [matplotlib](https://matplotlib.org/stable/users/getting_started/)
    * [pebble](https://github.com/noxdafox/pebble). The conda install can be found [here](https://anaconda.org/conda-forge/pebble)
    * [progressbar2](https://github.com/WoLpH/python-progressbar). The conda install can be found [here](https://anaconda.org/conda-forge/progressbar2).
    * [scipy](https://scipy.org/download/). The conda install can be found [here](https://anaconda.org/anaconda/scipy).
    * [pandas](https://pandas.pydata.org/docs/getting_started/install.html)
    * [biopython](http://biopython.org/wiki/Main_Page) 1.80 or later, see [note](#NOTE).
    * [psutil](https://github.com/giampaolo/psutil). The conda install can be found [here](https://anaconda.org/conda-forge/psutil). 
* [Exonerate](http://www.ebi.ac.uk/~guy/exonerate/) 2.40 or later
* [DIAMOND](https://github.com/bbuchfink/diamond/wiki). The conda install can be found [here](https://anaconda.org/bioconda/diamond).
* [BBtools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/). The conda install can be found [here](https://anaconda.org/bioconda/bbmap).



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

- HybPiper can now check target files for sequences with problematic low-complexity regions via command `hybpiper check_targetfile`.
- The `hybpiper assemble` command checks that the provided target file is formatted correctly and can be translated as expected (in the case of a nucleotide target file). Any issues are printed to screen and details are logged to file.
- The program DIAMOND can be used in place of BLASTX when mapping reads to target sequences.
- When using BLAST or DIAMOND, the `hybpiper stats` command now calculates the enrichment efficiency; previously this was only calculated when using BWA.
- The `hybpiper assemble` commands can now accept read files in compressed gzip format (suffix `*.gz`).
- Logging when running `hybpiper assemble` has been unified and extended to provide additional debugging information. A single log file is written per-sample in the sample directory e.g. `EG30_hybpiper_assemble_2021-12-02-10_45_56.log`. 
- Checks for all dependencies are now run by default when `hybpiper assemble` is run.
- Gene assemblies can no be performed using the SPAdes MDA (single cell)  mode.
- All Exonerate searches are now performed with the option `--refine full`; in the case of failure, a fallback run without this parameter is performed.
- In some situations HybPiper creates a gene sequence by 'stitching' together Exonerate hits from different SPAdes contigs. In some scenarios this can result in hit from different paralogs being joined together. HybPiper now performs a test (work in progress, currently sensitivity is low) to search for such 'chimeric' sequences, and provides warnings in the file `{sample_name}_genes_derived_from_putative_chimera_stitched_contig.csv` within each sample folder. Note that this chimera test is only performed in cases where a stitched contig has been created from multiple contigs, and paired-end reads are provided.
- By default, the SPAdes assembly folder is now deleted for each gene after contigs have been recovered. The user no longer needs to run `cleanup.py` after each run, and this script has been removed. Deleting the SPAdes directory dramatically reduces the total number of files produced by a completed run of HybPiper, which can be very useful when running it on and HPC with file number limits. To retain the SPAdes directory (i.e. for debugging purposes), the flag `--keep_intermediate_files` can be used.
- When running Intronerate (via flag `-run_intronerate` with command `hybpiper assemble`, a block of 10 'N' characters is inserted into supercontigs at locations where different SPAdes contigs have been concatenated. This behaviour can be turned off via the flag `--no_padding_supercontigs`.
- In cases where HybPiper recovers sequence for multiple non-contiguous segments of a gene, the gaps between the  segments will be padded with a number of 'N' characters. The number of Ns corresponds to the number of amino acids in 'best' protein reference for that gene that do not have corresponding SPAdes contig hits, multiplied by 3 to convert to nucleotides.
- The command `hybpiper stats` now writes pipeline run statistics directly to file, rather than to standard out.
- The command `hybpiper recover_sequences` now supports sequence recovery from a single sample.
- The command `hybpiper recover_sequences` now supports filtering of samples based on statistics generated via the `hybpiper stats` command (e.g. number of genes recovered, number of genes with paralogs, etc).
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
- `--bbmap_memory`. The amount of memory (RAM) in MB to use for BBmap when searching for chimeric stitched contigs. Default is 1000.
- `--chimeric_stitched_contig_edit_distance`. Minimum number of differences between one read of a read pair vs a stitched contig reference for a read pair to be flagged as discordant. Default is 5. 
- `--chimeric_stitched_contig_discordant_reads_cutoff`. Minimum number of discordant reads pairs required to flag a stitched contig as a potential chimera of contigs from multiple paralogs. Default is 5.
- `--keep_intermediate_files`. 1) If used, the SPAdes assembly folder for each gene will not be deleted. Note that previous versions of HybPiper retained the SPAdes assembly folders by default; they could previously be removed after running `reads_first.py` using the `cleanup.py` script. 2) If used, retain the `exonerate_hits.py` module log within each gene folder. Default behaviour is to delete each log after it has been copied to the main log file in the sample directory.
- `--distribute_hi_mem`. When used, distributing and writing reads to individual gene directories will be 40-50 percent faster, but can use more memory/RAM with large input files.
- `--single_cell_assembly`. Run SPAdes assemblies using MDA (single-cell) mode. 
- `--verbose_logging`. When used, the pipline logging will be much more verbose (particularly in the Exonerate stage), which can increase log file size dramatically.
- 

For subcommand **`hybpiper retrieve_sequences`**:

- `--single_sample_name`. Allows the user to specify a single sample name; sequences from this sample only will  be recovered in single `*.fasta` file.
- `--skip_chimeric_genes`. If used, skip retrieval of any genes/sample for which the stitiched contig was flagged as chimeric (potentially derived from multiple paralogs).
- `--fasta_dir`. User can specify a directory for output sequences.
- `--stats_file`. User can supply the stats fil produced by the command `hybpiper stats` for filtering samples (see below).
- `--filter_by`. User can filter samples based on statistics in the stats file.

For subcommand **`hybpiper stats`**:

- A positional parameter for the `sequence type` has been added, with the options {gene, supercontig}. Note that all length statistics are now provided as number of nucleotides, regardless of whether an amino-acid or nucleotide bait file was used.
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

- The parameter `--baitfile` or `-b` has been replaced with `--targetfile_dna/targetfile_aa` or `-t_dna/t-aa`, to more accurately describe the contents of the file provided. 
- The parameter `--length_pct` has been removed as a parameter to `hybpiper assemble` and is no longer used in internal code.
- The parameter `--thresh` (percent identity threshold for retaining Exonerate hits) for `hybpiper assemble` now defaults to 55 (previously 65).
- The flag `--check-depend` has been replaced by `hybpiper check_dependencies`. Dependency checking is now also performed every time the `hybpiper assemble` command is run.
- The flags `--no-blast`, `--no-distribute`, `--no-exonerate`, and `--no-assemble` have been replaced with the parameter `--start_from`. 
- The parameter `--timeout` has been replaced with `--timeout_assemble` and `--timeout_exonerate_contigs`, providing more fine-grained control over different stages of the assemble pipeline. 

***The following output files have been changed or removed***:

- ToDo

***The following output files/folders have been added***:

- The command `hybpiper paralog_retreiver` now write paralogs to two folders - one with all sequences, and the other without putative chimeric sequences.
- The command `hybpiper paralog_retreiver` now produces a heatmap showing paralogs detected.
- The command `hybpiper paralog_retreiver` now produces a report file in *tsv format listing genes with paralogs in above user-provided thresholds.
- The file containing paralog warnings (produced when multiple long contigs are present for a gene) has been renamed from `genes_with_paralog_warnings.txt` to `<sample_name>_genes_with_long_paralog_warnings.txt`. 
- In addition to the standard paralog warning produced when multiple long contigs are present for a gene, the `hybpiper assemble` command now provides a paralog warning when multiple short contigs are present which together cover the reference sequence for a given gene at a depth >1, across a given percentage length (default 75%) of the reference. These warning are written to each sample directory to the file `<sample_name>_genes_with_paralog_warnings_by_contig_depth.csv`
- `chimera_test_diagnostic_reads.sam`.
- XXX Other intronerate files, if we end up keeping them?
- The `hybpiper retrieve_sequences command` can now produce an output folder to contain fasta files.

    
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
