
import os,sys
from Bio import SeqIO

curdir = os.getcwd()

protein_directory = "/Users/mjohnson/Desktop/hybseq_pipeline/targets/mitochondrion"
alignment_directory = "/Users/mjohnson/mtdna"
outfilename = "mt_summary_mosses.txt"
#os.chdir(alignment_directory)
species = [d.split(".")[0] for d in os.listdir(alignment_directory) if os.path.isdir(os.path.join(alignment_directory, d)) and  not d.startswith(".")]
#outfile = open(outfilename,'w')
#os.chdir(protein_directory)
protein_list = ["%s.FAA" % i.split(".")[0] for i in os.listdir(protein_directory)]
#print sorted(protein_list)
sys.stdout.write("%s\t%s\n" % ("species","\t".join(a.split(".")[0] for a in protein_list)))
for sp in species:
	species_proteindir = os.path.join(alignment_directory,sp,"alignments/FAA")
	os.chdir(species_proteindir)
	sp_list = [len(SeqIO.read(protfile,'fasta').seq) if os.path.isfile(protfile) else 0 for protfile in protein_list]
	sys.stdout.write("%s\t%s\n" %
		(sp,"\t".join(str(a) for a in sp_list))
		)
#outfile.close()