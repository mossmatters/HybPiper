
import os,sys
from Bio import SeqIO

curdir = os.getcwd()

baitfilename = sys.argv[1]
sequence_directory = os.path.abspath(sys.argv[2])


species = [d.split(".")[0] for d in os.listdir(sequence_directory) if os.path.isdir(os.path.join(sequence_directory, d)) and  not d.startswith(".")]

protein_list = list(set([record.id.split("-")[1] + ".FAA" for record in SeqIO.parse(baitfilename,'fasta')]))
protein_list.sort()

#protein_list = ["%s.FAA" % i.split(".")[0] for i in os.listdir(protein_directory)]
sys.stdout.write("%s\t%s\n" % ("species","\t".join(a.split(".")[0] for a in protein_list)))
for sp in species:
    species_proteindir = os.path.join(sequence_directory,sp,"sequences/FAA")
    os.chdir(species_proteindir)
    sp_list = [len(SeqIO.read(protfile,'fasta').seq) if os.path.isfile(protfile) else 0 for protfile in protein_list]
    sys.stdout.write("%s\t%s\n" %
        (sp,"\t".join(str(a) for a in sp_list))
        )
#outfile.close()