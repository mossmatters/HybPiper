# HybPiper
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/hybpiper/README.html)

Current version: 2.1.8 (June 2024). See the change_log.md [here](https://github.com/mossmatters/HybPiper/blob/master/change_log.md)

[![DOI](https://zenodo.org/badge/6513/mossmatters/HybPiper.svg)](https://zenodo.org/badge/latestdoi/6513/mossmatters/HybPiper)

--
[Read our article in Applications in Plant Sciences (Open Access)](http://www.bioone.org/doi/full/10.3732/apps.1600016)

**HybPiper version 1.x** by Matt Johnson and Norm Wickett, Chicago Botanic Garden

If you would like to use HybPiper version 1.3, please download the `Final HybPiper 1.3 Version` [release](https://github.com/mossmatters/HybPiper/releases/tag/v1.3.1_final). The legacy wiki for version 1.3 can be found [here](https://github.com/mossmatters/HybPiper/wiki/HybPiper-Legacy-Wiki). 

**HybPiper version 2.0** by Matt Johnson (Texas Tech University) and Chris Jackson (Royal Botanic Gardens Victoria, Melbourne)


![hybpiper_logo](https://user-images.githubusercontent.com/55370301/168408947-e99a39e6-0d95-419a-8419-61b97315d863.png)

(logo by Elliot Gardner)


### Purpose

HybPiper was designed for targeted sequence capture, in which DNA sequencing libraries are enriched for gene regions of interest, especially for phylogenetics. HybPiper is a suite of Python scripts/modules that wrap and connect bioinformatics tools in order to extract target sequences from high-throughput DNA sequencing reads. 


Targeted bait capture is a technique for sequencing many loci simultaneously based on bait sequences. The HybPiper pipeline starts with high-throughput sequencing reads (for example from Illumina MiSeq), and assigns them to target genes using BLASTx/DIAMOND or BWA. The reads are distributed to separate directories, where they are assembled separately using SPAdes. The main output is a FASTA file of the (in frame) CDS portion of the sample for each target region, and a separate file with the translated protein sequence.

HybPiper also includes commands to extract the intronic regions flanking each exon, and investigate putative paralogs. For more information, [please see our wiki](https://github.com/mossmatters/HybPiper/wiki/).

HybPiper is run separately for each sample (single or paired-end sequence reads, with an optional file of unpaired reads in the latter scenerio). When HybPiper generates sequence 
files from the reads, it does so in a standardized directory hierarchy. Many of the post-processing commands rely on this directory hierarchy, so do not modify it after running the initial pipeline. It is a good idea to run the pipeline for each sample from the same directory. You will end up with one directory per run of HybPiper, and some of the later commands take advantage of this predictable directory structure.


---
# Dependencies
* [Python](https://www.python.org/downloads/) 3.7 or later (see [note](#NOTE)), along with the Python libraries:
    * [seaborn](https://seaborn.pydata.org/installing.html)
    * [matplotlib](https://matplotlib.org/stable/users/getting_started/)
    * [pebble](https://github.com/noxdafox/pebble). The conda install can be found [here](https://anaconda.org/conda-forge/pebble)
    * [progressbar2](https://github.com/WoLpH/python-progressbar). The conda install can be found [here](https://anaconda.org/conda-forge/progressbar2).
    * [scipy](https://scipy.org/download/). The conda install can be found [here](https://anaconda.org/anaconda/scipy).
    * [pandas](https://pandas.pydata.org/docs/getting_started/install.html)
    * [biopython](http://biopython.org/wiki/Main_Page) 1.80 or later, see [note](#NOTE).
    * [psutil](https://github.com/giampaolo/psutil). The conda install can be found [here](https://anaconda.org/conda-forge/psutil). 
* [Exonerate](http://www.ebi.ac.uk/~guy/exonerate/) 2.40 or later
* [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)  2.9.0 
* [DIAMOND](https://github.com/bbuchfink/diamond/wiki). The conda install can be found [here](https://anaconda.org/bioconda/diamond).
* [BWA](http://bio-bwa.sourceforge.net/)
* [BBtools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/). The conda install can be found [here](https://anaconda.org/bioconda/bbmap).
* [SPAdes](http://bioinf.spbau.ru/en/spades) 3.15.0
* [GNU Parallel](http://www.gnu.org/software/parallel/)
* [samtools](http://www.htslib.org/) 1.14
* [MAFFT](https://mafft.cbrc.jp/alignment/software/). ***New requirement for HybPiper version 2.1.0***

#### NOTE: 

***Biopython***

Biopython 1.80 is required as it contains bug fixes for the SearchIO module that HybPiper 2.0 uses to parse output from Exonerate.

**Update 01/12/2022:** Biopython 1.80 has now been released, and is available as a conda install [here](https://anaconda.org/conda-forge/biopython).

Prior to this release, HybPiper users needed to install an in-progress version by cloning the repository from https://github.com/biopython/biopython and following the instructions under the heading 'Installation from Source' [here](https://biopython.org/wiki/Download). Alternatively, the HybPiper conda package for all versions prior to 2.1.0 includes pre-release Biopython version 1.80.dev0; note that this **causes problems** if the user then installs Biopython 1.79 via conda in the HybPiper conda environment. The HybPiper conda package for version 2.1.0 now uses the conda Biopython 1.80 release package. 

***Python and SPAdes***

The SPAdes assembler prior to version 3.15.4 is incompatible with Python 3.10 and 3.11, see [here](https://github.com/ablab/spades/issues/873). The highest SPAdes version [available](https://anaconda.org/bioconda/spades) as a macOS conda package is 3.15.2. This means that the macOS conda package for HybPiper is restricted to Python version 3.9 at the highest. This restriction will be removed once SPAdes >= 3.15.4 is released as a macOS conda package.  

---
# Setup

We strongly recommend installing HybPiper using [conda](https://docs.conda.io/en/latest/miniconda.html) with a new environment. This will install HybPiper, all required Python packages, and all required external programs. If you have conda installed and the channels `bioconda` and `conda-forge` have already been added, this can be done using the command:

```
conda create -n hybpiper hybpiper
```

...followed by:

```
conda activate hybpiper
```

For full installation instructions, including details on how to install on **Macs with Apple Silicon (M1/M2/M3 chips)**, please see our wiki page:

[https://github.com/mossmatters/HybPiper/wiki/Installation](https://github.com/mossmatters/HybPiper/wiki/Installation)


Once all dependencies are installed you can download a test dataset [here](https://github.com/mossmatters/HybPiper/raw/master/test_dataset.tar.gz), extract the `*.tar.gz` file, and execute the `run_hybpiper_test_dataset.sh` script from the `test_dataset` directory for a demonstration of HybPiper.


----

# Pipeline Input

Full instructions on running the pipeline, including a step-by-step tutorial using a small test dataset, are available on our wiki:

[https://github.com/mossmatters/HybPiper/wiki](https://github.com/mossmatters/HybPiper/wiki)

### High-Throughput DNA Sequencing Reads

Before running the pipeline, you will likely want "cleaned" FASTQ file(s) - one or two depending on whether your sequencing was single or paired-end. Reads should have adapter sequences removed and should be trimmed to remove low quality base calls.

### Target Sequences

You will also need to construct a "target" file of gene regions. The target file should contain one gene region per sequence, with exons "concatenated" into a contiguous sequence. For more information on constructing the target file, see the wiki, or view the example file in: `test_dataset/test_targets.fasta`

There can be more than one "source sequence" for each gene in the target file. This can be useful if the target enrichment baits were designed from multiple sources - for example a transcriptome in the focal taxon and a distantly related reference genome.

----

# Pipeline Output

HybPiper will map the reads to the target sequences, sort the reads by gene, assemble the reads for each gene separately, align the contigs to the target sequence, and extract a coding sequence from each gene. Output from each of these phases is saved in a standardized directory hierarchy, making it easy for post-processing commands to summarize information across many samples.

For example, the coding sequence for gene "gene001" for sample "EG30" is saved in a FASTA file:

`EG30/gene001/EG30/sequences/FNA/gene001.FNA`

...and a list of genes for which a sequence could be extracted can be found:

`EG30/genes_with_seqs.txt`

For a full description of HybPiper output, [see the wiki](https://github.com/mossmatters/HybPiper/wiki).


-----
# Changelog

**2.1.8** *June, 2024*

Added the new subcommand `hybpiper filter_by_length`. Used to filter the sequence output of `hybpiper retrieve sequences` by absolute length and/or length relative to mean length in target file representatives. This is done on a per-sample/per-gene basis, rather than the sample-level filtering available in `hybpiper retrieve_sequences`. See [wiki](https://github.com/mossmatters/HybPiper/wiki#hybpiper-filter_by_length) for more information. For a full list of changes see the [changelog](https://github.com/mossmatters/HybPiper/blob/master/change_log.md).

**2.1.0** *December, 2022*

Added the new subcommand `fix_targetfile`. Assists with filtering target file sequences based on length and sequence complexity, and in ensuring sequences are in the correct reading frame (the latter for nucleotide target files only). For a full list of changes see the [changelog](https://github.com/mossmatters/HybPiper/blob/master/change_log.md).

***New dependencies***:

- MAFFT


**2.0 release candidate** *May, 2022*

This update involves a substantial refactor of the HybPiper pipeline, with changes to the internal code, additional functionality, and additional output. For a full list of changes see the [changelog](https://github.com/mossmatters/HybPiper/blob/master/change_log.md).

***New dependencies***: 

- Python 3.6 or later, along with the Python libraries:
   - seaborn
   - matplotlib
   - pebble
   - progressbar2
   - scipy
   - pandas
   - BioPython 1.80
   - psutil
- DIAMOND
- BBtools (BBmap.sh, BBmerge.sh)

    
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

---

# Citation
Johnson, M. G., Gardner, E. M., Liu, Y., Medina, R., Goffinet, B., Shaw, A. J., Zerega, N. J. C, and  Wickett, N. J. (2016). HybPiper: Extracting Coding Sequence and Introns for Phylogenetics from High-Throughput Sequencing Reads Using Target Enrichment. Applications in Plant Sciences, 4(7), 1600016. doi:10.3732/apps.1600016

