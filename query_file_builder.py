#!/usr/bin/env python

# Assumes that the directory contains files of Hyb-Seq bait nucleotide sequences.
# 
# For each file:
#     Use exonerate to determine the best hit to the genome assembly.
#     Output only that protein to standard out.
    
import sys, os, subprocess,argparse,shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def blastx(proteinfilename,assemblyfilename,evalue,cleanup,prefix):
    """Conducts a blastx search of the assembly against a database of all the proteins
    Selects the best hit for each gene, and returns a list of proteins to save by id.
    Assumes there is a preformatted blastx database in the same directory as the proteinfilename."""
    
    blastx_cmd = "blastx -query %s -db %s -evalue %f -outfmt 6 > %s/temp.blastx" %(assemblyfilename,proteinfilename,evalue,prefix)
    #print "Blasting, please wait"
    proc = subprocess.call(blastx_cmd, shell=True)
    
    #print "Reading BLASTx results."
    protein_hash = {}
    blastx_results = open("%s/temp.blastx"%prefix)
    for line in blastx_results.readlines():
        line = line.split()
        protid = line[1].split("-")
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
    if cleanup:os.remove("%s/temp.blastx"%prefix)
    return protein_hash
                    
def makeblastdb(proteinfilename,cleanup):
    if os.path.exists(proteinfilename+".psq"):
        return
    else:
        makeblastdb_cmd = "makeblastdb -in %s -dbtype prot >makeblastdb.log" % proteinfilename
        proc = subprocess.call(makeblastdb_cmd,shell=True)
    
    if cleanup: os.remove("makeblastdb.log")
    return                    

    
def main():
    parser = argparse.ArgumentParser(description="query_file_builder.py; Generate 'Tailored' bait files for each protein. REQUIRES: Blast Command Line Tools")
    parser.add_argument("proteinfile",help="""FASTA file containing several copies of each protein. The sequence IDs should be of the format: 
        >speciesname-proteinname""")
    parser.add_argument("assemblyfile",help="""FASTA file of DNA sequence assembly.""")
    parser.add_argument("--evalue",help="e-value cutoff for BLAST search. default = 0.001",default=0.001,type=float) 
    parser.add_argument("--full_cleanup",help="Delete all temporary files. default = no",action="store_true")
    parser.add_argument("--overwrite",help="Overwrite the bait file if it already exists. default = no",action="store_true")
    
    args=parser.parse_args()    
    proteinfilename = args.proteinfile
    shutil.copy(proteinfilename,".")
    proteinfilename =  os.path.basename(proteinfilename)
    assemblyfilename = args.assemblyfile
    
    prefix = os.path.basename(assemblyfilename).split(".")[0]
    if not os.path.exists(prefix):
        os.makedirs(prefix)
    shutil.copy(assemblyfilename,"%s/assembly.fasta"%prefix)
    
    outfilename = "%s/baitfile.FAA"%prefix
    if os.path.exists(outfilename):
        if not args.overwrite:
            print("Bait file %s already exists!"  % outfilename)
            return

    outfile = open(outfilename,'w')
    
    makeblastdb(proteinfilename,args.full_cleanup)
    besthits = blastx(proteinfilename,assemblyfilename,args.evalue,args.full_cleanup,prefix)
    all_proteins = SeqIO.to_dict(SeqIO.parse(proteinfilename,'fasta'))
    for prot in besthits:
        seqID = "%s-%s" % (besthits[prot][0],prot)
        SeqIO.write(all_proteins[seqID],outfile,'fasta')
    outfile.close()
    


if __name__ == "__main__": main()    