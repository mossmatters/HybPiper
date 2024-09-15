"""AIM: The goal of this piece of code is to test the performance of HybPiper
as the test sequences are more and more different from the target sequences.
The motivation is to explore scenarios where the species we are using is 
very diverged from the ones used to generate the baits of the HybSeq array.

This pipeline uses files within test_dataset directory in HybPiper. We will create 
a loop and each time we are going to randomly introduce variation in the test_sequences.

Every time we will generate heatmap and summary stats to explore how HybPiper performance changes.

Notes for me:
We need to be extracareful to keep the HybPiper hierarchical directory structure, otherwise
it will fail.
Every new file should be created within the test_dataset folder.

Remember to change R code accordingly  when we do the LOOP!
Is mutable module or class?
"""
## PENDING: Generate the loop across iterations with appropriate variables

# We first want to run HybPiper's pipeline as it is
Rscript multiSampleLoop.R

# Now we generate heatmap and summstats files.
python ../get_seq_lengths.py test_targets.fasta namelist.txt dna > test_seq_lengths.txt
Rscript gene_recovery_heatmap.R
python hybpiper_stats.py test_seq_lengths.txt namelist.txt > ~/HybPiper/PAFTOL_pipeline/test_stats.txt

# Now we change the test reads fastq. We plan to use the mutable module or class.
python ~/HybPiper/PAFTOL_pipeline/mutateFastq.py 

