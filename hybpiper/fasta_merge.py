#!/usr/bin/env python

helptext ='''This script will take a list of FASTA files and concatenate them for use in 
phylogenetic inference. The sequence headers (up until the first space) must be identical
in each individual FASTA file.

Individual gene sequences should be aligned prior to running this script!

This script requires BioPython to read/write FASTA sequences.'''

import os,sys,argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def read_sequences(fastafiles):
    '''Given a list of FASTA file names, read in each sequence to a dictionary of dictionaries, one per file'''
    return {filename:SeqIO.to_dict(SeqIO.parse(filename,'fasta')) for filename in fastafiles}

def get_unique_names(gene_dict):
    '''Given the dictionary of SeqRecord dictionaries, return a list of the unique sequence headers'''
    all_names = []
    for gene in gene_dict:
        all_names += list(gene_dict[gene].keys())
    return set(all_names)

def insert_sequences(gene_dict,unique_names):
    '''Given the dictionary of dictionaries, insert blank sequences if any are missing for a gene'''
    inserted_sequences = 0
    for gene in gene_dict:
        for name in unique_names:
            if name not in gene_dict[gene]:
                gene_length = len(next(iter(gene_dict[gene].values())))
                gene_dict[gene][name] = SeqRecord(Seq("-"*gene_length),id=name)
                inserted_sequences += 1
    sys.stderr.write("{} Empty sequences inserted across all genes.\n".format(inserted_sequences))            
    return gene_dict

def concatenate_sequences(gene_dict,fastafiles,unique_names):
    '''Given a dictionary of dictionaries with complete sampling in each gene, write out concatenated sequences to stdout. Returns a list of partition lengths.'''    
    new_seq_dict = {}
    partition_lengths = []
    for gene in fastafiles:
        for name in unique_names:
            try:
                new_seq_dict[name] += gene_dict[gene][name]
            except KeyError:
                new_seq_dict[name] = gene_dict[gene][name]
        partition_lengths.append(len(next(iter(gene_dict[gene].values()))))
    for final_seq in new_seq_dict:
        SeqIO.write(new_seq_dict[final_seq],sys.stdout,'fasta')            
    final_seq_length = len(new_seq_dict[final_seq])
    sys.stderr.write("Final conatenated sequence length: {}\n".format(final_seq_length))
    return partition_lengths

def raxml_partition(fastafiles,partition_lengths,partition_type):
    '''Generate a raxml partition file for the given fastafiles. User specifies the partition type'''
    gene_start = 1
    partition_file = open("partition.raxml",'w')
    
    if partition_type == 'CODON':
        for g in range(len(fastafiles)):
            codon3_start = gene_start + 2
            codon3_end = gene_start + partition_lengths[g] - 1
            codon1_end = codon3_end - 2
            codon2_start = gene_start + 1
            codon2_end = codon3_end - 1
            partition_file.write("{},{}{}={}-{}\\3,{}-{}\\3\n".format("DNA",fastafiles[g],"12",gene_start,codon1_end,codon2_start,codon2_end))
            partition_file.write("{},{}{}={}-{}\\3\n".format("DNA",fastafiles[g],"3",codon3_start,codon3_end))
            gene_start = codon3_end + 1
    else:
        for g in range(len(fastafiles)):
            gene_end = gene_start + partition_lengths[g] - 1
            partition_file.write("{},{}={}-{}\n".format(partition_type,fastafiles[g],gene_start,gene_end))
            gene_start = gene_end + 1
        partition_file.close()    




def main():
    parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--fastafiles",nargs='+',help="List of Fasta Files. Can use wildcard on Linux/Mac systems")
    parser.add_argument("--filelist",help="File containing list of Fasta files. Alternative to --fastalist")
    parser.add_argument("--raxml",help="Create a partition file 'partitions.raxml' intended for raxml in the current directory. For amino acid sequences, select the substitution model. To specify a separate model for 1st/2nd vs. 3rd codon positions, select CODON.",
        choices = ['DNA','WAG','JTT','CODON'
                    ],default=None)
        
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    
    if args.fastafiles:
        #print args.fastafiles
        if args.filelist:
            sys.stderr.write("Specify either a list of FASTA files or a file containing names, not both!\n")
            sys.exit(1)
        else:
            fastafiles = args.fastafiles    
        
    elif args.filelist:
        #print args.filelist
        if os.path.isfile(args.filelist):
            fastafiles = [x.rstrip() for x in open(args.filelist)]
        else:
            sys.stderr.write("File containing list of FASTA files not found!")
            sys.exit(1)    
    
    else:
         sys.stderr.write("You must specify the FASTA files as a list or in a file.\n")
         sys.exit(1)

    sys.stderr.write("{} FASTA files found.\n".format(len(fastafiles)))
    gene_dict = read_sequences(fastafiles)
    
    sys.stderr.write("All sequences read successfully.\n")    
    unique_names = get_unique_names(gene_dict)
    sys.stderr.write("{} Unique names found. If you were expecting fewer sequences, check your IDs!\n".format(len(unique_names)))
    gaps_inserted = insert_sequences(gene_dict,unique_names)

    partition_lengths = concatenate_sequences(gaps_inserted,fastafiles,unique_names)
    
    if args.raxml:
        raxml_partition(fastafiles,partition_lengths,args.raxml)

if __name__ == "__main__":main()