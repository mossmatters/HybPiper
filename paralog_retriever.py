#!/usr/bin/env python

helptext = '''This script will retrieve paralog nucleotide (CDS) sequences for a specified
gene in all samples located in namelist.txt. It writes all the (unaligned) sequences to stdout.
If a sample does not have paralogs for that gene, the sequence in the FNA directory is retrieved instead.'''

import os,sys,argparse
from Bio import SeqIO

def retrieve_seqs(path,name,gene):
    seqs_to_write = None
    if os.path.isdir(os.path.join(path,name,gene,name,'paralogs')):
        seqs_to_write = [x for x in SeqIO.parse(os.path.join(path,name,gene,name,'paralogs','{}_paralogs.fasta'.format(gene)),'fasta')]
        num_seqs = str(len(seqs_to_write))
    elif os.path.isfile(os.path.join(path,name,gene,name,'sequences','FNA','{}.FNA'.format(gene))):
        seqs_to_write = SeqIO.read(os.path.join(path,name,gene,name,'sequences','FNA','{}.FNA'.format(gene)),'fasta')
        num_seqs = "1"
    
    if seqs_to_write:
        SeqIO.write(seqs_to_write,sys.stdout,'fasta')
    else:
        num_seqs = "0"
    
    return num_seqs
            

def main():
    parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('namelist',help='Text file containing list of HybPiper output directories, one per line.')
    parser.add_argument('gene',help="Name of gene to extract paralogs")
    
    args = parser.parse_args()
    
    namelist = [x.rstrip() for x in open(args.namelist)]
    num_seqs = []
    for name in namelist:
        path, name = os.path.split(name)
        if not name:
            path,name = os.path.split(path)
        num_seqs.append(retrieve_seqs(path,name,args.gene))
    sys.stderr.write("{}\t{}\n".format(args.gene,"\t".join(num_seqs)))
    
if __name__ == "__main__":main()
