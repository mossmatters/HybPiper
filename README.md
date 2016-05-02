#HybPiper
--
*Manuscript in Revision*

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




