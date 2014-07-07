#!/usr/bin/env python

import sys,os,errno,argparse
from Bio import SeqIO

helptext = """

usage: python distribute_targets.py baitfile\n


Given a file containing all of the "baits" for a target enrichment, create separate
FASTA files with all copies of that bait. Multiple copies of the same bait can be 
specified using a "-" delimiter. For example, the following will be sorted in to the same
file:

Anomodon-rbcl
Physcomitrella-rbcl

Output directories can also be created, one for each target category
	(the default is to put them all in the current one)
The field delimiter may also be changed.

NOTE: This script will APPEND sequences to existing files. 
 If you are re-running this script, delete the existing output first!
"""


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
        
def distribute_targets(baitfile,dirs,delim):
	targets = SeqIO.parse(baitfile,'fasta')
	
	for prot in targets:
		#Get the 'basename' of the protein
		prot_cat = prot.id.split(delim)[-1]
		
		if dirs:
			mkdir_p(prot_cat)        
		
		outfile = open(os.path.join(prot_cat,"{}_baits.fasta".format(prot_cat)),'a')	
		if prot.id.startswith("Physco"):
			SeqIO.write(prot,outfile,'fasta')
		outfile.close()

def help():
	print helptext
	return
		
def main():
	parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("--no_dirs",help="Do not generate separate directories for each protein-- output all to the current directory.", action="store_true",default=False)
	parser.add_argument("-d","--delimiter",help="Field separating FASTA ids for multiple sequences per target. Default is '-' . For no delemeter, write None", default="-")
	parser.add_argument("baitfile",help="FASTA file containing bait sequences")
	args = parser.parse_args()
	
	if args.no_dirs:
		distribute_targets(args.baitfile,dirs=None,delim=args.delimiter)
	else:
		distribute_targets(args.baitfile,dirs=True,delim=args.delimiter)
	


if __name__=="__main__":main()