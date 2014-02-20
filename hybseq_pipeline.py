#!/usr/env python

# Given an assembly file of genomic DNA and a file containing target proteins:
# 	1. Use exonerate to determine hits for each contig.
# 	2. Load contigs into a Biopython seqRecord dictionary.
# 	3. For each protein hit, create a FASTA file with all contigs.
# 		a. The contigs will be in order of their hit location to the target proteins.

######REQUIREMENTS########
# Python 2.6 or later
# exonerate in your $PATH
# Biopython
##########################

import sys, os, subprocess,math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

id_threshold = 55 #Percent identity between protein and contig hit.
first_search_filename = "exonerate_results.fasta"

def initial_exonerate(proteinfilename, assemblyfilename):
	"""Conduct exonerate search, returns a dictionary of results.
	Using the ryo option in exonerate, the header should contain all the useful information."""
	
	outputfilename = "exonerate_results.fasta"
	exonerate_ryo = '">%ti,%qi,%qab,%qae,%pi,(%tS),%tab,%tae\\n%tcs\\n"'
	exonerate_command = "exonerate -m protein2genome --showalignment no --showvulgar no -V 0 --ryo %s %s %s >%s" % (exonerate_ryo,proteinfilename,assemblyfilename,outputfilename)
	print exonerate_command
	#print exonerate_ryo
	#proc = subprocess.Popen(['exonerate','-m','protein2genome','--showalignment','no','-V','0','--showvulgar','no','--ryo',exonerate_ryo,proteinfilename,assemblyfilename])
	proc = subprocess.call(exonerate_command,shell=True)
	protHitsCount = 0
	#proc.wait()
	records = SeqIO.to_dict(SeqIO.parse("exonerate_results.fasta",'fasta'))
	#proc.stdout.close()
	
	return records

def protein_sort(records):
	"""Given the Biopython dictionary, return a dictionary of proteins indexed by their hits.""" 
	proteinHits = {}
	for contig in records:
		hit = records[contig].id.split(",")
		protein = hit[1]
		if protein in proteinHits:
			proteinHits[protein]["assemblyHits"].append(",".join(hit))
			proteinHits[protein]["hit_start"].append(int(hit[2]))
			proteinHits[protein]["hit_end"].append(int(hit[3]))
			proteinHits[protein]["percentid"].append(float(hit[4]))
			proteinHits[protein]["hit_strand"].append(hit[5][1])
			proteinHits[protein]["target_begin"].append(int(hit[6]))
			proteinHits[protein]["target_end"].append(int(hit[7]))
		else:
			proteinHits[protein] = {"assemblyHits" : [",".join(hit)],
					"hit_start" : [int(hit[2])],
					"hit_end" : [int(hit[3])],
					"percentid" : [float(hit[4])],
					"hit_strand" : [hit[5][1]],
					"target_begin" : [int(hit[6])],
					"target_end" : [int(hit[7])],
					"name" : protein
					}
	return proteinHits

def get_contig_order(prot):
	"""Given the dictionary of hits for a protein, return the dictionary with the fields sorted by start location."""
	sorting_order = sorted(range(len(prot["hit_start"])),key=lambda k:prot["hit_start"][k])
	
	prot["assemblyHits"] = [prot["assemblyHits"][i] for i in sorting_order]
	prot["hit_start"] = [prot["hit_start"][i] for i in sorting_order]
	prot["hit_end"] = [prot["hit_end"][i] for i in sorting_order]
	prot["percentid"] = [prot["percentid"][i] for i in sorting_order]
	prot["hit_strand"] = [prot["hit_strand"][i] for i in sorting_order]
		
	return prot

def filter_by_percentid(prot,thresh):
	"""Given a protein dictionary, return a protein dictionary minus entries with percentID below a threshold"""
	kept_indicies = [i for i in range(len(prot["percentid"])) if prot["percentid"][i] > thresh]
	return keep_indicies(kept_indicies,prot)		

def supercontig_exonerate(supercontig,protseq):
	"""Given a long, joined contig and a protein sequence, return the exonerate hit(s)"""
	exonerate_ryo = '>%ti,%qi,%qab,%qae,%pi,(%tS)\\n%tcs\\n'
	SeqIO.write(protseq,"temp.prot.fa",'fasta')
	SeqIO.write(supercontig,"temp.contig.fa",'fasta')
	print "Conducting exonerate search"
	proc = subprocess.Popen(['exonerate','-m','protein2genome','--showalignment','no','-V','0','--showvulgar','no','--ryo',exonerate_ryo,"temp.prot.fa","temp.contig.fa"],stdout=subprocess.PIPE)
	proc.wait()
	#print proc.stdout.readlines()
	supercontig_cds = SeqIO.parse(proc.stdout,'fasta')
	
	#print [i.id for i in supercontig_cds]
	return supercontig_cds

def sort_byhitloc(seqrecord):
	"""Key function for sorting based on the start location of a hit record."""
	return int(seqrecord.id.split(",")[2])
	
def fullContigs(prot,sequence_dict,assembly_dict,protein_dict):
	"""Generates a contig from all hits to a protein. 
	If more than one hit, conduct a second exonerate search with the original contigs
	stitched together."""
	numHits = len(prot["assemblyHits"])
	sequence_list = []
	contigHits = []
	
	#print numHits
	if numHits == 1:
		print prot["assemblyHits"]
		print prot["percentid"]
		return str(sequence_dict[prot["assemblyHits"][0]].seq)	#If only one hit to this protein.
	else:
		for hit in range(len(prot["assemblyHits"])):
			assembly_seq_name = prot["assemblyHits"][hit].split(",")[0]
			
			if assembly_seq_name not in contigHits:		#Only add each contig once.
				if prot["hit_strand"][hit] == "+":
					sequence_list.append(assembly_dict[assembly_seq_name])
				else:
					sequence_list.append(assembly_dict[assembly_seq_name].reverse_complement())
			contigHits.append(assembly_seq_name)
	print prot["name"]
 	print [i for i in prot["assemblyHits"]]
 	print [(prot["hit_start"][i],prot["hit_end"][i]) for i in range(len(prot["hit_start"]))]
 	print prot["hit_strand"]
 	print prot["percentid"]
	supercontig = SeqRecord(Seq("".join(str(b.seq) for b in sequence_list)),id=prot["name"])
	#print supercontig
	
	#Need to remove contigs if they have the same basename
	
	supercontig_cds = supercontig_exonerate(supercontig,protein_dict[prot["name"]])
	#Sort the supercontigs by hit location to the protein.
	joined_supercontig_cds = [b for b in supercontig_cds]
	joined_supercontig_cds.sort(key=sort_byhitloc)
	print joined_supercontig_cds
	return str(Seq("".join(str(b.seq) for b in joined_supercontig_cds)))	
	#print joined_supercontig_cds
	#print ""
	return joined_supercontig_cds


def find_longest_hit(prot):
	"""Given a protein dictionary, determine the assembly hit with the longest sequence"""
	max_hit_length = 0
	max_hit_loc = 0
	for i in range(len(prot["hit_start"])):
		hit_length = abs(int(prot["hit_start"][i]) - int(prot["hit_end"][i]))
		if hit_length > max_hit_length:
			hit_length = max_hit_length
			max_hit_loc = i
	return max_hit_loc

def keep_indicies(kept_indicies, prot):
	"""Given a list of indicies to keep and a protein dictionary, return the dictionary with only the specified entries remaining"""
	assHit = []
	hitstart = []
	hitend = []
	percentid = []
	strands = []
	targetbegin =[]
	targetend =[]

	for a in kept_indicies:
		assHit.append(prot["assemblyHits"][a])
		hitstart.append(prot["hit_start"][a])
		hitend.append(prot["hit_end"][a])
		percentid.append(prot["percentid"][a])
		strands.append(prot["hit_strand"][a])
		targetbegin.append(prot["target_begin"][a])
		targetend.append(prot["target_end"][a])

	prot["assemblyHits"] = assHit
	prot["hit_start"] = hitstart
	prot["hit_end"] = hitend
	prot["percentid"] = percentid
	prot["hit_strand"] = strands
	prot["target_begin"] = targetbegin
	prot["target_end"] = targetend

	return prot

def overlapping_contigs(prot):
	"""Given a protein dictionary, determine whether the hit ranges are overlapping,
	and save only those contigs that are not completely subsumed by other contigs."""
	range_list = [(prot["hit_start"][i],prot["hit_end"][i]) for i in range(len(prot["hit_start"]))]
	kept_indicies = range_connectivity(range_list)
	return keep_indicies(kept_indicies,prot)

def range_connectivity(range_list):
	"""Given two sorted lists, representing the beginning and end of a range,
	Determine "connectivity" between consecutive elements of the list.
	For each connected segment, determine whether one segement "subsumes" the other."""
	starts = [a[0] for a in range_list]
	ends = [a[1] for a in range_list]

	connected_ranges = []
	collapsed_ranges = []
	num_breaks = 0
	for i in range(len(starts)):
		try:
			next_start = starts[i+1]
		except IndexError:	#Special case at the end of the list
			if starts[i] <= ends[i-1]:
				connected_ranges.append((starts[i],ends[i]))
				collapsed_ranges.append(subsume(connected_ranges))
			else:
				num_breaks += 1
				connected_ranges.append((starts[i],ends[i]))
				collapsed_ranges.append(subsume(connected_ranges))
				connected_ranges = []
			continue
		if next_start <= ends[i]:
			connected_ranges.append((starts[i],ends[i]))
		else:
			num_breaks += 1
			connected_ranges.append((starts[i],ends[i]))
			collapsed_ranges.append(subsume(connected_ranges))
			connected_ranges = []
	if False: #num_breaks == 0:
		kept_indicies = [range_list.index(i) for i in connected_ranges]
		return kept_indicies
	else:
		#List contains other lists, need to flatten this to just tuples.
		flattened_list = []
		for a in range(len(collapsed_ranges)):
			if isinstance(collapsed_ranges[a], list):
				for i in collapsed_ranges[a]:
					flattened_list.append(i)
			else:
				flattened_list.append(collapsed_ranges[a])
		kept_indicies = [range_list.index(i) for i in flattened_list]
		return kept_indicies
	
def subsume(connected_ranges):
	"""Determine whether one range in the list of ranges completely encompasses the other."""
	starts = [a[0] for a in connected_ranges]
	ends = [a[1] for a in connected_ranges]
	for i in range(len(connected_ranges)):
		if starts[i] == min(starts) and ends[i] == max(ends):
			return connected_ranges[i]
	return connected_ranges
def tuple_overlap(a,b):
	"""Given two tuples of length two, determine if the ranges overlap"""
	return a[0] < b[0] < a[1] or b[0] < a[0] < b[1]

def reciprocal_best_hit(prot,proteinHits):
	"""Given a protein dictionary and the dictionary of all protein dictionaries,
		Return the protein dictionary minus any contigs that have higher percentage hits to other proteins."""
	protname = prot["name"]
	kept_indicies=[]
	for contig in prot["assemblyHits"]:
		contigname = contig.split(",")[0]
		contig_idx = prot["assemblyHits"].index(contig)
		maxProt = protname
		for otherProt in proteinHits:
			#print "checking %s vs %s" %(protname, proteinHits[otherProt]["name"])
			otherprot_contiglist = [x.split(",")[0] for x in proteinHits[otherProt]["assemblyHits"]]
			if proteinHits[otherProt]["name"] != protname:
				if contigname in otherprot_contiglist:
					full_contigname = [b for b in proteinHits[otherProt]["assemblyHits"] if contigname in b][0]
					print contig, full_contigname
					otherHit_idx = proteinHits[otherProt]["assemblyHits"].index(full_contigname)
					
					target_ranges = [sorted((prot["target_begin"][contig_idx],prot["target_end"][contig_idx])),sorted((proteinHits[otherProt]["target_begin"][otherHit_idx],proteinHits[otherProt]["target_end"][otherHit_idx]))]
					print target_ranges
					#Check that the two contig hits have overlapping ranges.
					if tuple_overlap(target_ranges[0],target_ranges[1]): 				
						print prot["percentid"][contig_idx],proteinHits[otherProt]["percentid"][otherHit_idx]
						if prot["percentid"][contig_idx] < proteinHits[otherProt]["percentid"][otherHit_idx]:
							print "contig %s is a better hit to %s" %(contigname,otherProt)
							maxProt = proteinHits[otherProt]["name"]
					else:
						print "ranges did not overlap"
		if maxProt == protname:
			kept_indicies.append(contig_idx)
	return keep_indicies(kept_indicies,prot)		



def myTranslate(nucl):
	"""Given a raw sequence of nucleotides, return raw sequence of amino acids."""
	#print nucl
	nucseq = Seq(nucl)
	#print nucseq
	aminoseq = nucseq.translate()
	return str(aminoseq)

def help():
	print "USAGE: python hybseq_pipeline.py proteinfile assemblyfile prefix"
	print "The program Exonerate must be in your $PATH."
	print "You must have BioPython installed"
	print "A protein and a nucleotide directory will be created in the current directory with the prefix."
	return	

def main(): 
	if len(sys.argv) < 4:
		help()
		return
		
	proteinfilename = sys.argv[1]
	assemblyfilename = sys.argv[2]
	prefix = sys.argv[3]

	try:
		proteinfile = open(proteinfilename)
	except IOError:
		print "The file %s could not be opened!" %proteinfilename
		return()
		
	try:
		assemblyfile = open(assemblyfilename)
	except IOError:
		print "The file %s could not be opened!" % assemblyfilename
		return()
	assembly_dict = SeqIO.to_dict(SeqIO.parse(assemblyfile,'fasta'))
	protein_dict = SeqIO.to_dict(SeqIO.parse(proteinfile,'fasta'))
	
	if os.path.exists(first_search_filename): 	#Shortcut for Testing purposes
		print "Reading exonerate results from file."
		sequence_dict = SeqIO.to_dict(SeqIO.parse(first_search_filename,'fasta'))
	else:
		print "Starting exonerate search, please wait."
		sequence_dict = initial_exonerate(proteinfilename,assemblyfilename)
	proteinHits = protein_sort(sequence_dict)

	print "There were %i exonerate hits." %	len(sequence_dict)
	print "There were %i unique proteins hit." % len(proteinHits)
	
	directory_name = "%s/FNA" % prefix
	if not os.path.exists(directory_name):
		os.makedirs(directory_name)

	directory_name = "%s/FAA" % prefix
	if not os.path.exists(directory_name):
		os.makedirs(directory_name)
	
 	for prot in proteinHits:
 		#Put contigs in order along the protein.
 		print proteinHits[prot]["name"]
 		print "Initial hits: ", proteinHits[prot]["assemblyHits"]
 		proteinHits[prot] = get_contig_order(proteinHits[prot])
 		print "After get_contig_order: ", proteinHits[prot]["assemblyHits"]
		#Remove contigs that are suboptimal hits. Only one protein hit allowed per contig.
		proteinHits[prot] = reciprocal_best_hit(proteinHits[prot],proteinHits)
 		print "After RBH: " , proteinHits[prot]["assemblyHits"]
 		if len(proteinHits[prot]["assemblyHits"]) == 0:
 			continue		#All hits have been filtered out
 		
 		#Filter out contigs with a hit below a threshold
 		proteinHits[prot] = filter_by_percentid(proteinHits[prot],id_threshold)
 		print "After filter_by_percent_id:", proteinHits[prot]["assemblyHits"]
 		if len(proteinHits[prot]["assemblyHits"]) == 0:
 			continue		#All hits have been filtered out
 		
 		#Delete contigs if their range is completely subsumed by another hit's range.
 		proteinHits[prot] = overlapping_contigs(proteinHits[prot])
 		print "After overlapping_contigs: " , proteinHits[prot]["assemblyHits"]
 		#Stitch together a "supercontig" containing all the hits and conduct a second exonerate search.	
 		if len(proteinHits[prot]["assemblyHits"]) == 0:
 			continue		#All hits have been filtered out

 		nucl_sequence = fullContigs(proteinHits[prot],sequence_dict,assembly_dict,protein_dict)
 		#print nucl_sequence
		print proteinHits[prot]["name"]
		#print nucl_sequence
		amino_sequence = myTranslate(nucl_sequence)

		amino_filename = "%s/FAA/%s.FAA" % (prefix,prot)#.split("-")[1])
		amino_file = open(amino_filename,'w')
		amino_file.write(">%s\n%s\n" % (prefix,amino_sequence))
		amino_file.close()
		
		nucleo_filename = "%s/FNA/%s.FNA" % (prefix,prot)#.split("-")[1])
		nucleo_file = open(nucleo_filename,'w')
		nucleo_file.write(">%s\n%s\n" % (prefix,nucl_sequence))
		nucleo_file.close()
# 			
	
	proteinfile.close()
	assemblyfile.close()


if __name__ == "__main__":main()