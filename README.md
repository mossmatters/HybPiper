#HybPiper

[![Join the chat at https://gitter.im/mossmatters/HybPiper](https://badges.gitter.im/mossmatters/HybPiper.svg)](https://gitter.im/mossmatters/HybPiper?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)  [![DOI](https://zenodo.org/badge/6513/mossmatters/HybPiper.svg)](https://zenodo.org/badge/latestdoi/6513/mossmatters/HybPiper)

--
[Read our article in Applications in Plant Sciences (Open Access)](http://www.bioone.org/doi/full/10.3732/apps.1600016)

by Matt Johnson and Norm Wickett, Chicago Botanic Garden

![](examples/hybpiper_logo.png)

(logo by Elliot Gardner)


###Purpose

HybPiper was designed for targeted sequence capture, in which DNA sequencing libraries are enriched for gene regions of interest, especially for phylogenetics. HybPiper is a suite of Python scripts that wrap and connect bioinformatics tools in order to extract target sequences from high-throughput DNA sequencing reads. 


Targeted bait capture is a technique for sequencing many loci simultaneously based on bait sequences.
HybPiper pipeline starts with high-throughput sequencing reads (for example from Illumina MiSeq), and assigns them to target genes using BLASTx or BWA.
The reads are distributed to separate directories, where they are assembled separately using SPAdes. 
The main output is a FASTA file of the (in frame) CDS portion of the sample for each target region, and a separate file with the translated protein sequence.

HybPiper also includes post-processing scripts, run after the main pipeline, to also extract the intronic regions flanking each exon, investigate putative paralogs, and calculate sequencing depth. For more information, [please see our wiki](https://github.com/mossmatters/HybPiper/wiki/).

HybPiper is run separately for each sample (single or paired-end sequence reads). When HybPiper generates sequence files from the reads, it does so in a standardized directory hierarchy. Many of the post-processing scripts rely on this directory hierarchy, so do not modify it after running the initial pipeline. It is a good idea to run the pipeline for each sample from the same directory. You will end up with one directory per run of HybPiper, and some of the later scripts take advantage of this predictable directory structure.


---
#Dependencies
* Python 2.7 or later (to use the argparse module for help documents)
* [BIOPYTHON 1.59 or later](http://biopython.org/wiki/Main_Page) (For parsing and handling FASTA and FASTQ files)
* [EXONERATE](http://www.ebi.ac.uk/~guy/exonerate/) (For aligning recovered sequences to target proteins)
* [BLAST command line tools](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (Aligning reads to target protiens)
* [SPAdes](http://bioinf.spbau.ru/en/spades) (Assembling reads into contigs)
* [GNU Parallel](http://www.gnu.org/software/parallel/) (Handles parallelization of searching, assembling, and aligning)

*Required for BWA version of the pipeline and for the intron and depth calculation scripts*:

* [BWA](http://bio-bwa.sourceforge.net/) (Aligns reads to target nucleotide sequences)
* [samtools](http://www.htslib.org/) (Read/Write BAM files to save space).

**NOTE:** A previous version of the pipeline required Velvet and CAP3 for assembly. These have been unreliable at assembling individual genes, and SPAdes has replaced them.

---
#Setup
We have successfully installed HybPiper on MacOSX and Linux (Centos 6). All of the bioinformatics tools can be installed with [homebrew](brew.sh) or [linuxbrew](linuxbrew.sh).

For full installation instructions, please see our wiki page:

[https://github.com/mossmatters/HybPiper/wiki/Installation](https://github.com/mossmatters/HybPiper/wiki/Installation)

Once all dependencies are installed, execute the `run_tests.sh` script from the `test_dataset` directory for a demonstration of HybPiper.


----

#Pipeline Input

Full instructions on running the pipeline, including a step-by-step tutorial using a small test dataset, is available on our wiki:

[https://github.com/mossmatters/HybPiper/wiki](https://github.com/mossmatters/HybPiper/wiki)

###High-Throughput DNA Sequencing Reads

Before running the pipeline, you will need "cleaned" FASTQ file(s)-- one or two depending on whether your sequencing was single or paired-end. Reads should have adapter sequences removed and should be trimmed to remove low quality base calls.

###Target Sequences

You will also need to construct a "target" file of gene regions. The target file should contain one gene region per sequence, with exons "concatenated" into a contiguous sequence. For more information on constructing the target file, see the wiki, or view the example file in: `test_dataset/test_targets.fasta`

There can be more than one "source sequence" for each gene in the target file. This can be useful if the target enrichment baits were designed from multiple sources-- for example a transcriptome in the focal taxon and a distantly related reference genome.

----

#Pipeline Output

HybPiper will map the reads to the target sequences, sort the reads by gene, assemble the reads for each gene separately, align the contigs to the target sequence, and extract a coding sequence from each gene. Output from each of these phases is saved in a standardized directory hierarchy, making it easy for post-processing scripts to summarize information across many samples.

For example, the coding sequence for gene "gene001" for sample "EG30" is saved in a FASTA file:

`EG30/gene001/EG30/sequences/FNA/gene001.FNA`

and a list of genes for which a sequence could be extracted can be found:

`EG30/genes_with_seqs.txt`

For a full description of HybPiper output, [see the wiki](https://github.com/mossmatters/HybPiper/wiki).


-----
#Changelog

**1.2** *January, 2017*

* Features 

	- Added `--unpaired` flag. When using paired-end sequencing reads, a third read file may be specified with this flag. Reads will be mapped to targets separately, but will be used along with paired reads in contig assembly.

	- Added `--target` flag. Using this flag, HybPiper will only use targets with the specified name in the Alignment/Exon Extraction phase. All other targets for that locus will only be used in the Mapping/Read Sorting phase.

	- Added `--timeout` flag, which uses GNU Parallel to kill processes (i.e. Spades or Exonerate) if they take X percent longer than average. Use if there are a lot of stuck jobs (`--timeout 1000`)

* Bug Fixes

 - Can accommodate Solexa FASTQ paired headers
 - Fixed `spades_runner.py` not recognizing `--cpu` on redos
 - Prints more meaningful messages for some common errors
 - Can accommodate `prefix` not being in current directory
 - Deletes sorted reads on restart to prevent double counting reads.


**1.1** *May, 2016*: Release associated with manuscript in *Applications in Plant Sciences*.

-- Added paralog_investigator.py, which detects and extracts long exons from putative paralogs in all genes in one sample.

-- Added paralog_retriever.py, which retrieves sequences generated by paralog_investigator.py for many samples (or the coding sequence generated by exonerate_hits.py if no paralogs are detected).

-- Added a test_dataset of 13 genes for 9 samples, and a shell script to run the test data through the main script and several post-processing scripts.

-- Fixed bug involving calling HybPiper with a relative path such as: `../reads_first.py`

-- `reads_first.py --check_depend` now checks for SPAdes, BWA, and Samtools

-- Full revision of README, which is now shorter. Full tutorials on installing and running HybPiper are now on the Wiki.



**1.0** *Feb, 2016*: Initial fully-featured release associated with submission of manuscript to *Applications in Plant Sciences*.

-- Sequence assembly now uses SPAdes rather than Velvet + CAP3

---

#Citation
Johnson, M. G., Gardner, E. M., Liu, Y., Medina, R., Goffinet, B., Shaw, A. J., Zerega, N. J. C, and  Wickett, N. J. (2016). HybPiper: Extracting Coding Sequence and Introns for Phylogenetics from High-Throughput Sequencing Reads Using Target Enrichment. Applications in Plant Sciences, 4(7), 1600016. doi:10.3732/apps.1600016

