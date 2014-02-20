#!/usr/env python

# Assumes that the directory contains files of Hyb-Seq bait nucleotide sequences.
# 
# For each file:
# 	Use exonerate to determine the best hit to the genome assembly.
# 	Output only that protein to standard out.
	
import sys, os, subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def exonerate(proteinfilename,assemblyfilename):
	"""Return the best exonerate result for a protein."""
	exonerate_ryo = '%qi\t%ti\t%s\n'
	proc = subprocess.Popen(['exonerate','-m','cdna2genome','--showalignment','no','-V','0','--showvulgar','no','-n','1','--percent','50','--ryo',exonerate_ryo,proteinfilename,assemblyfilename],stdout=subprocess.PIPE)
	proc.wait()
	results = [r.split() for r in proc.stdout.readlines()]
	scores = [int(r[2]) for r in results]
	print scores
	if len(scores) > 1:
		best_result = results[scores.index(max(scores))][0]	
		return best_result
	else:
		return None

def usearch(proteinfilename,assemblyfilename):
	usearch_cmd = "usearch -usearch_global %s -db %s -id 0.8 -strand both -blast6out temp.usearch -quiet" %(proteinfilename,assemblyfilename)
	proc = subprocess.Popen(usearch_cmd,shell=True)
	filter_cmd = "cat temp.usearch | sort -n -t$'\t' -k3 |tail -1|cut -f1 > temp.name"
	proc = subprocess.Popen(filter_cmd,shell=True)
	
	
	best_result = open("temp.name").readline()
	
	
	return best_result
	
def blastx(proteinfilename,assemblyfilename):
	"""Conducts a blastx search of the assembly against a database of all the proteins
	Selects the best hit for each gene, and returns a list of proteins to save by id.
	Assumes there is a preformatted blastx database in the same directory as the proteinfilename."""
	
	blastx_cmd = "blastx -query %s -db %s -evalue 0.001 -outfmt 6 > temp.blastx" %(assemblyfilename,proteinfilename)
	#print "Blasting, please wait"
	proc = subprocess.call(blastx_cmd, shell=True)
	
	#print "Reading BLASTx results."
	protein_hash = {}
	blastx_results = open("temp.blastx")
	for line in blastx_results.readlines():
		line = line.split()
		protid = line[1].split(",")
		tax = protid[0]
		protname = protid[1]
		score = line[2]
		if protname in protein_hash:
			if score > protein_hash[protname][1]:
				protein_hash[protname] = (tax,score)
		else:
			protein_hash[protname] = (tax,score)
	#print protein_hash
	blastx_results.close()
	return protein_hash
					
				

	
def main():
	#bait_directory = sys.argv[1]
	proteinfilename = sys.argv[1]
	assemblyfilename = sys.argv[2]
	#curdir = os.getcwd()
	#os.chdir(bait_directory)
	#protein_list = os.listdir(bait_directory)
	
	besthits = blastx(proteinfilename,assemblyfilename)
	all_proteins = SeqIO.to_dict(SeqIO.parse(proteinfilename,'fasta'))
	for prot in besthits:
		seqID = "%s-%s" % (besthits[prot][0],prot)
		SeqIO.write(all_proteins[seqID],sys.stdout,'fasta')
		
		
# 	for prot in protein_list:
# 		#besthit = exonerate(prot,assemblyfilename)
# 		besthit = usearch(prot,assemblyfilename)
# 		print prot.split(".")[0], besthit
# 		#os.remove("temp.usearch")
# 		#os.remove("temp.name")


if __name__ == "__main__": main()	