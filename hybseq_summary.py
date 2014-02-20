
import os
from Bio import SeqIO

curdir = os.getcwd()

protein_directory = "/scratch/mjohnson/hyseqPipeline/targets/plastid"
alignment_directory = "/scratch/mjohnson/hyseqPipeline/assemblies/plastid/new_plastid/pipeline"
outfilename = "/scratch/mjohnson/hyseqPipeline/cp_summary_mosses.txt"
os.chdir(alignment_directory)
species = os.listdir(alignment_directory)

outfile = open(outfilename,'w')
os.chdir(protein_directory)
protein_list = [i.split(".")[0] + ".FAA" for i in os.listdir(protein_directory)]
print sorted(protein_list)
outfile.write("%s\t%s\n" % ("species","\t".join(a.split(".")[0] for a in protein_list)))
for sp in species:
	species_proteindir = os.path.join(alignment_directory,sp,"FAA")
	os.chdir(species_proteindir)
	sp_list = [len(SeqIO.read(protfile,'fasta').seq) if os.path.isfile(protfile) else 0 for protfile in protein_list]
	outfile.write("%s\t%s\n" %
		(sp,"\t".join(str(a) for a in sp_list))
		)
outfile.close()