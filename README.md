#HybPiper

*Manuscript in Preparation*

by Matt Johnson and Norm Wickett, Chicago Botanic Garden

![](examples/hybpiper_logo.png)

*Purpose* 

Targeted bait capture is a technique for sequencing many loci simultaneously based on bait sequences.
This pipeline starts with Illumina reads, and assigns them to target genes using BLASTx or BWA.
The reads are distributed to separate directories, where they are assembled separately using SPAdes. 
The main output is a FASTA file of the (in frame) CDS portion of the sample for each target region, and a separate file with the translated protein sequence.

An optional script run after the main pipeline can also extract the 
intronic regions flanking each exon.

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
##Install the dependencies.
Normal installations are fine for most tools, as long as they are in your `$PATH`.

For Python and BioPython, the [Anaconda distribution](https://www.continuum.io/downloads) may be easiest to use and install. 

You can check the installation by executing the `readsfirst.py` script with the `--check-depend` flag.

For Spades, if you download pre-compiled binaries, make sure that the contents of `spades/bin` and `spades/share/spades` are in your `$PATH`

###MACOSX installation notes

Both Exonerate and Velvet require zlib. 
Perhaps the easiest way to install Exonerate is with homebrew: http://brew.sh/

Install homebrew, and then "tap" the science repository:

	brew tap homebrew/science

Now install exonerate:

	brew install exonerate

Homebrew will install zlib as part of the dependencies for exonerate. 

You can also install samtools, spades, and bwa this way:

	brew install bwa
	brew install samtools
	brew install spades

For velvet, this command line worked for me in Mac OS 10.9:

`make 'MAXKMERLENGTH=67' LDFLAGS='-L third-party/zlib-1.2.3/ -lm’`

##Preparing your files
**IMPORTANT**: If you are using BLAST to map reads to targets, you need a PROTEIN baitfile.

If you are using BWA, you need a NUCLEOTIDE baitfile. 

Construct a "bait file" of protein sequences in FASTA format. It is ok to have more than one "source sequence" per bait. 
The ID for each sequence should include the bait source and the protein ID, separated by a hyphen. For example:

```
>Amborella-atpH
MNPLISAASVIAAGLAVGLASIGPGVGQGTAAGQAVEGIARQPEAEGKIRGTLLLSLAFM
>Aneura-atpH
MNPLIPAASVIAAGLAVGLASIGPGIGQGTAAGQAVEGIARQPEAEGKIRGTLLSSPASM
>Amborella-rbcL
MSPKTETKASAGFKAGVKDYRLTYYTPDYETLATDILAAFRVTPQPGVPPEEAGAAVAAE
>Aneura-rbcL
MSPQTETKAGVGFKAGVKDYRLTYYTPEYETKETDILAAFRMTPQPGVPPEEAGAAVAAE
```

##Running the pipeline

HybPiper is run separately for each sample (single or paired-end sequence reads). When HybPiper generates sequence files from the reads, it does so in a standardized directory heirarchy. Many of the post-processing scripts rely on this directory heirarchy, so do not modify it after running the initial pipeline. It is a good idea to run the pipeline for each sample from the same directory. You will end up with one directory per run of HybPiper, and some of the later scripts take advantage of this predictable directory structure.

To execute the entire pipeline, create a directory containing the paired-end read files.
The script `reads_first.py` will create a directory based on the fastq filenames (or use the `--prefix` flag):

`Anomodon-rostratus_L0001_R1.fastq` ---> `Anomodon-rostratus/`


The following command will execute the entire pipeline on a pair of Illumina read files, using the baits in the file `baits.fasta`. The HybPiper scripts should be in a different directory:

```python /Users/mehmattski/HybPiper/reads_first.py -r MySpecies_R1.fastq MySpecies_R2.fastq -b baits.fasta```

The BLASTx version of the pipeline (default) will map the reads to amino acid bait sequences sequences.
Although it is slower than the BWA version, it may have higher specificity. Reads may not align to divergent nucleotide bait sequences, which are required for the BWA version.
If you find the recovery efficiency is poor with BWA, you may want to try again with BLASTx.


#Pipeline Scripts

###`reads_first.py`
A wrapper script that:

1. Can check if all dependencies are installed correctly. (`--check-depend`)

2. Creates sub-directories.

3. Calls all downstream analyses

You can tell the script to skip upstream steps (for example: `--no-blast`) but the script will still assume that the output files of these steps still exist!

This script will call, in order:

1. Blastx (or BWA)

2. `distribute_reads_to_targets.py` (or `distribute_reads_to_targets_bwa.py`)  and `distribute_targets.py`

3. Run SPAdes assembler.

4. `exonerate_hits.py`

Some program-specific options may be passed at the command line. 

For example, the e-value threshold for BLASTX (`--evalue`, default is `1e-9`), the coverage-cutoff level for SPAdes assemblies (`--cov_cutoff`, default is `8`), or the percent-identity threshold for aligning contigs to the reference sequences (`--thresh`, the default is `55`).

Use `reads_first.py -h` for a full list.	

###`distribute_reads_to_targets.py`
After a BLASTx search of the raw reads against the target sequences, the reads need to be 
sorted according to the successful hits. This script takes the BLASTx output (tabular)
and the raw read files, and distributes the reads into FASTA files ready for assembly.

If there are multiple BLAST results (for example, one for each read direction),
concatenate them prior to sorting.

###`distribute_reads_to_targets_bwa.py`
After a BWA search of the raw reads against the target sequences, the reads need to be 
sorted according to the successful hits. This script takes the bam output (parsed using samtools)
and the raw read files, and distributes the reads into FASTA files ready for assembly.

If there are multiple BAM results (for example, one for each read direction),
concatenate them prior to sorting.

###`distribute_targets.py`
Given a file containing all of the "baits" for a target enrichment, create separate
FASTA files with all copies of that bait. Multiple copies of the same bait can be 
specified using a "-" delimiter. For example, the hits against the following two baits will be sorted in to the same file:

`Anomodon-rbcl`

`Physcomitrella-rbcl`

Given multiple baits, the script will choose the most appropriate 'reference' sequence
using the highest cumulative BLAST scores across all hits. If the search was BWA rather than BLAST, it will use the bait with the highest total BWA alignment score for each gene.

Output directories can also be created, one for each target category (the default is to put them all in the current one). The field delimiter may also be changed.


###`exonerate_hits.py`
This script generates alignments of SPAdes contigs against the target sequence. 

If BLASTx is used, the model is `protein2genome`

If BWA is used, the model is `coding2genome`

Contigs are not expected to overlap much. An inital exonerate search is filtered for hits that are above a certain threshold (default is 55, can be changed in ``reads_first.py`` with `--thresh`). Contigs that pass this filter are arranged in order along the alignment. If one contig completely subsumes the range of another contig, the longer contig is saved. 

To maximize the chance that Exonerate may find introns, all contigs that pass the previous steps are concatenated into a "supercontig" and the exonerate search is repeated. Once again, unique hits that pass a percent identity threshold and do not overlap with longer hits are saved, and from this the full length CDS and protein sequences are generated.

In the older version of the script ("assemble-first"), the script was used after the `query_file_builder` is complete. The minimal inputs are the tailored bait file and the assembly.

If run immediately after query_file_builder, use the --prefix flag to specify the file names:

EXAMPLE COMMAND LINE:

```exonerate_hits.py speciesName/baitfile.FAA speciesName/assembly.fasta --prefix=speciesName```

The threshold for accepting an exonerate hit can be adjusted with `-t` (Default: 55 percent)

-----

#RESULTS AND OUTPUT FILES

##Base directory
(Specified by ``--prefix`` or generated from the read file names)
 
1. The master bait file is copied here

2. A BLAST or BWA database is generated

3. One directory is generated for every gene with BLAST or BWA hits. 

4. A file "bait_tallies.txt" summarizes which bait sources were chosen.
	
##Gene Directory
1.	SPAdes. Final assembly is at "GeneName_contigs.fasta"

2.	Fasta file for reference bait chosen by the `distribute_targets.py` script.

3.	Directory of Exonerate results (with same name as the sample)

##Exonerate Directory

1.	`exonerate_results.fasta` -- Results of the initial exonerate search for all contigs.

2.	`supercontig_exonerate.fasta` -- Long concatenated contig from final exonerate search.

3.	sequences directory

##Sequences directory

1.	`FNA/GeneName.FNA`: In-frame nucleotide sequence

2.	`FAA/GeneName.FAA`: Amino acid sequence.
	
#Summary

The major steps of the pipeline include:

1. Blast (or BWA) search of the reads against the target sequences.

2. Distribution of reads into separate directories, one per gene.

3. Assembly of reads for each gene into contigs with SPAdes.

4. Conduct one or more exonerate searches for each contig in the assembly. If multiple contigs match the same protein in non-overlapping sequences, stitch the hits together into a “supercontig”

5. In a subdirectory, generate separate FASTA files containing either the nucleotide (FNA) or amino acid (FAA) sequence for each protein. 

-----
#Paralogs

In many HybSeq bait designs, care is taken to avoid enrichment for genes with paralogous sequences in the target genomes. However, gene duplication and paleopolyploidy (especially in plants) makes it difficult to completely avoid paralogs. However, given enough read depth, it may be possible to identify paralogous sequences with HybPiper. 

If SPAdes assembler generates multiple contigs that contain coding sequences represeting 75% of the length of the reference protein, HybPiper will print a warning for that gene. It will also print the names of the contigs in ```prefix/genename/paralog_warning.txt``` and will print a list of all genes with paralog warnings to ```prefix/genes_with_paralog_warnings.txt```.

If many of your genes have paralogs, one approach could be to add the paralog coding sequence to your bait file as a separate gene, and re-running HybPiper. Reads that have a better mapping to the paralog will be sorted accordingly.

-----

#Introns

Frequently, probe sequences for HybSeq are designed from coding regions (exons) only. However, given the read lengths now attained by some sequencing technologies (e.g. Illumina MiSeq 2x300), there is the potential to recover flanking intron and intergenic regions. Introns are more likely to be variable with in species or species complexes.

The figure below shows the depth of MiSeq reads aligned to the draft genome of *Artocarpus altilis*. The gray lines represent depth within a sliding window (50 bp) across the genome scaffold for 22 *Artocarpus* samples. The dark line is the overall average depth, and the red bars represent the exons. Substantial depth is acheived up to 400 bp away from the exons for most samples.

We have extended HybPiper to extract sequences flanking the coding sequence for each gene.

###`intronerate.py`

Given a completed run of `reads_first.py` for a sample, run this script to generate "gene" sequences for each locus. The script will generat two new sequence files for each gene:

**supercontig**: A sequence containing all assembled contigs with a unique alignment to the reference protein, concatenated into one sequence.

**introns**: The supercontig with the exon sequences removed.

	python interonerate.py --prefix hybseq_directory
	
Specify the name of a directory generated by ```reads_first.py``` in the prefix argument.

The default behavior is to refer to the file `genes_with_seqs.txt` to recover full length sequences only for the genes where exons were previously recovered. You may optionally supply a file containing a list of genes with `--genelist filename`

**NOTE**: The script will extract all sequence *NOT* annotated as exons by exonerate. This may be introns (or intergenic sequence), but it may also be mis-assembled contigs. While it may be difficult ot tell whether the sequence is "real" from a single sample, I recommend running `intronerate.py` on several samples. Then, extract the supercontig sequences with `retrieve_sequences.py` and align them. Sequences that appear in only one sample are probably from mis-assembled contigs and may be trimmed, for example using Trimal.

----


#After the pipeline
Optional utilities after running the pipeline for multiple assemblies: 

**NOTE**: for these utilities to work, the files must be in the same directory hierarchy created by the pipeline. (i.e. `species/sequences/FAA/` and `species/sequences/FNA/`)

###`cleanup.py`

HybPiper generates a lot of output files. Most of these can be discarded after the run. This script handles deleting unnecessary files, and can reduce the size of the directory created by HybPiper by 75%.

####Example Command Line

```
python cleanup.py hyb_seq_directory
```

By default the script will delete all the files generated by velvet. Other options may be added in the future.


###`get_seq_lengths.py`

This script will summarize the recovery of genes from multiple samples. If you have all of these separate runs of the HybPiper in the same directory, create a `namelist.txt` file that contains a list of all the HybPiper directories for each sample (one per line):

```
Sample1
Sample2
Sample3
```

Specify the location of the bait file and whether it is amino acid or nucleotide on the command line:

####Example Command Line

`python get_seq_lengths.py baitfile.fasta namelist.txt dna > gene_lengths.txt`

The script will output a table to `stdout`. The first line is a header containing the gene names. The second line contains the length of the reference for each gene. If there are multiple reference sequences for each gene, an average is reported. The remaining lines are the lengths recovered by the HybPiper for each sample, one sample per line (one column per gene). If the gene is missing, a 0 is indicated.

A warning will print to stderr if any sequences are longer than 1.5x the average reference length for that gene.

####Example output
```
Species	26S	18S
MeanLength	3252.0	1821.6
funaria	3660	1758
timmia	3057	1821
```

###`gene_recovery_heatmap.R`

This script takes the ouput of `get_seq_lengths.py` and creates a figure to visualize the recovery efficiency.

Unlike the python scripts, you will need to open the R script in a text editor or RStudio and edit a few parameters before running the script within R. 

The script requires two R packages: `gplots` and `heatmap.plus` 
Install these using `install.packages` before running the script.
You will need to set the name of your input file (the one produced by `get_seq_lengths.py`) at the top of the script.
The output will look something like this:

![heatmap](examples/plastids_heatmap.png)

Each row shows a sample, and each column is a gene (in this case, one of the 44 chloroplast genes). The amount of shading in each box corresponds to the length of the gene recovered for that sample, relative to the length of the reference (bait). 

In this case, there are a few samples for which few or no genes were recovered (white rows) and a few genes that were not recovered in any sample (white columns). 



###`retrieve_sequences.py`

This script fetches the sequences recovered from the same gene for many samples and generates an unaligned multi-FASTA file for each gene. 

This script will get the sequences generated from multiple runs of the HybPiper (reads_first.py).
Have all of the runs in the same directory (sequence_dir). 
It retreives all the gene names from the bait file used in the run of the pipeline.

####Example Command Line

`python retrieve_sequences.py baitfile.fasta sequence_dir dna`

You must specify whether you want the protein (aa) or nucleotide (dna) sequences.

If you ran `intronerate.py` on these samples, you can also specify "supercontig" or "intron" to recover those sequences instead.

Will output unaligned fasta files, one per gene, to the current directory.



-----

#DEPRECATED SCRIPTS
These scripts are left over from a version of the pipeline that started with sequence assemblies, rather than raw reads.

##`query_file_builder.py`

This script generates the "tailored baitfile" for the species by choosing the best representative at each gene using a BLAST search against the assembly file. It also sets up all of the necessary file hierarchy to run the next step of the pipeline.

The input to the script requires a fasta file containing protein bait sequences, as described in Setup, and the nucleotide assembly file.

The script will use the prefix of the assembly file to generate a directory containing all the results. For the cleanest results, create a new directory, and use relative or absolute paths to indicate the locations of both the protein and assembly file.

For example, if one level up there is one directory containing baits and another containing assemblies:

EXAMPLE COMMAND LINE

`query_file_builder.py ../baits/all_plastid_baits.FAA ../assemblies/speciesName.fasta`



