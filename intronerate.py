#!/usr/bin/env python

helptext = '''
This script will take the output of a run of HybSeqPipeline (exon sequences) and attempt
to extract intron sequences from the velvet assemblies. It is important that the
directory structure of HybSeqPipeline was not disturbed, so that it can be used to collect
information and re-execute exonerate.
'''

import sys,os,argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from exonerate_hits import range_connectivity,tuple_overlap


def get_contig_info():
    '''Given the prefix of a run of HybSeqPipeline, retreive the stats from exonerate'''
    statspath = "exonerate_stats.csv"
    contig_info = [x.rstrip().split(',') for x in open(statspath).readlines()]
    #print contig_info
    sorted_contig_stats = sorted(contig_info, key = lambda x: int(x[2]))
    #print sorted_contig_stats
    return sorted_contig_stats
    
def make_intron_supercontig(contig_info,gene,prefix):
    cap3contigs = SeqIO.to_dict(SeqIO.parse("../{}_contigs.fasta".format(gene),'fasta'))
    intron_supercontig = SeqRecord(Seq(''))
    for i in contig_info:
        if i[5] == "(+)":
            intron_supercontig += cap3contigs[i[0]]
        elif i[5] == "(-)":
            intron_supercontig += cap3contigs[i[0]].reverse_complement()    
        else:
            sys.stderr.write("Strandedness not found!")
            sys.exit(1)    
    intron_supercontig.id = '{}-{}'.format(prefix,gene)
    intron_supercontig.description = ''
    SeqIO.write(intron_supercontig,'sequences/intron/{}_supercontig.fasta'.format(gene),'fasta')    
    
def re_run_exonerate(gene,target=new_faa):
    if target == new_faa:
        exonerate_cmd = "exonerate -m protein2genome -q sequences/FAA/{}.FAA -t sequences/intron/{}_supercontig.fasta --verbose 0 --showalignment no --showvulgar no --showtargetgff yes > intronerate_raw.gff".format(gene, gene)
    else:
        exonerate_cmd = "exonerate -m protein2genome -q ../{}_baits.fasta -t sequences/intron/{}_supercontig.fasta --verbose 0 --showalignment no --showvulgar no --showtargetgff yes > intronerate_raw.gff".format(gene, gene)
    sys.stderr.write("[CMD] {}\n".format(exonerate_cmd))
    os.system(exonerate_cmd)

def parse_gff(filename):
    '''Parse a GFF file created for a single gene, return a list of lists containing the annotation info'''
    with open(filename) as gff_file:
        gff_dump = gff_file.read()
        gff_split = gff_dump.split("# --- END OF GFF DUMP ---")
        raw_hits = [x.split('\n') for x in gff_split[:-1]]
        hits = []
        for h in raw_hits:
            new_hit = []
            for line in h:
                if not line.startswith("#") and not line == "":
                    new_hit.append(line.rstrip().split('\t'))
            hits.append(new_hit)
    return hits


def filter_gff(hits):
    hits_to_keep = []
    hits = sorted(hits,key= lambda x: int(x[0][3]))
    #Get only the features annotated as genes
    gene_annotations = [x for y in hits for x in y if x[2] == 'gene']
    #Get the start,end, and score for each gene annotation
    range_list = [(int(x[3]),int(x[4])) for x in gene_annotations]
    kept_indicies = range_connectivity(range_list)
    if len(kept_indicies) > 1:
        overlapping_indicies = []
        non_overlapping_indicies = []
        for ix in range(len(kept_indicies)-1):
            #print range_list[ix],range_list[ix+1],tuple_overlap(range_list[ix],range_list[ix+1])
            if tuple_overlap(range_list[ix],range_list[ix+1]):
                overlapping_indicies.append((kept_indicies[ix],kept_indicies[ix+1]))
            else:
                non_overlapping_indicies.append(kept_indicies[ix])
        for pair in overlapping_indicies:
            if int(gene_annotations[pair[0]][5]) > int(gene_annotations[pair[1]][5]):
                non_overlapping_indicies.append(pair[0])
            else:
                non_overlapping_indicies.append(pair[1])
        if not tuple_overlap(range_list[-2],range_list[-1]):
                non_overlapping_indicies.append(kept_indicies[-1])
        return     [hits[x] for x in non_overlapping_indicies]    
                
    else:
        return [hits[x] for x in kept_indicies]            
    #print [gene_annotations[x] for x in kept_indicies]

def get_new_gff(kept_hits):
    flatter_list = []
    for hit in kept_hits:
        for line in hit:
            flatter_list.append("\t".join(line))
    return "\n".join(flatter_list) + '\n'    

def remove_exons(gff_filename,supercontig_filename,mode="all"):
    '''Given a supercontig and corresponding annotation, remove the exon sequences. In "intron" mode, only return sequences specifically annotated as introns'''
    exon_starts = []
    exon_ends = []
    gff = open(gff_filename).readlines()
    for line in gff:
        line = line.rstrip().split("\t")
        if len(line) > 2:
            if line[2] == "exon":
                exon_starts.append(int(line[3]))
                exon_ends.append(int(line[4]))
    supercontig = SeqIO.read(supercontig_filename,'fasta')
    exonless_contig = SeqRecord(Seq(''),id=supercontig.id)
    start = 0
    for exon in range(len(exon_starts)):
        exonless_contig += supercontig[start:exon_starts[exon]-1] 
        start = exon_ends[exon]
    exonless_contig += supercontig[start:]    
    exonless_contig.description = ''
    return exonless_contig    

def check_for_files(gene,prefix):
    '''Check to see if the files needed for intronerate are really present'''
    if os.path.isfile("{}/{}/exonerate_stats.csv".format(gene,prefix)):
        if os.path.isfile("{}/{}/sequences/FAA/{}.FAA".format(gene,prefix,gene)):
            if os.path.isfile("{}/{}_contigs.fasta".format(gene,gene)):
                return True
            else:
                return False
        else:
            return False
    else:
        return False
def main():
    parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--genelist",help="Optional list of genes to retreive coverage. Default is to use genes_with_seqs.txt")
    parser.add_argument("--prefix",help="Prefix of sample directory generated by HybSeqPipeline",required=True)
    parser.add_argument("--no-exonerate",help = "Don't re-run exonerate, use existing intronerate gff files.",action='store_true',default=False)
    parser.add_argument("--use_target",help="Align the supercontig to the original target sequences, rather than the newly generated FAA",default=False,action="store_true")
    #parser.add_argument("--introns-only",help = "In the intron.fasta file for each gene, only write regions annotated as introns by exonerate. Default: all non-exon regions are written to introns.fasta.",action="store_true",default=False)
    args=parser.parse_args()
    
    if len(sys.argv) < 2:
        print(helptext)
        sys.exit(1)

    if os.path.isdir(args.prefix):
        os.chdir(args.prefix)
        basedir = os.getcwd()
        prefix = os.path.split(basedir)[1]
    else:
        sys.stderr.write("Directory {} not found!\n".format(args.prefix))    

    if args.genelist:
        genelist = [x.split()[0] for x in open(args.genelist).readlines()]
    else:
        genelist = [x.split()[0] for x in open('genes_with_seqs.txt').readlines()]

    with open("intron_stats.txt",'w') as intron_stats_file:    
        full_gff = ''
        for gene in genelist:
            if check_for_files(gene,prefix):
                os.chdir("{}/{}".format(gene,prefix))
                contig_info = get_contig_info()
                if not os.path.exists("sequences/intron"):
                    os.makedirs("sequences/intron")
                make_intron_supercontig(contig_info,gene,prefix)
                if not args.no_exonerate:
                    if args.use_target:
                        re_run_exonerate(gene,target=bait)
                    else:
                        re_run_exonerate(gene,target=new_faa)
                hits = parse_gff("intronerate_raw.gff")
                kept_hits = filter_gff(hits)
                with open("intronerate.gff",'w') as new_gff:
                    new_gff_string = get_new_gff(kept_hits)
                    num_introns = new_gff_string.count("intron\t")
                    intron_stats_file.write("{}\t{}\t{}\n".format(prefix,gene,num_introns))
                    sys.stderr.write("{} introns found for {}.\n".format(num_introns,gene))
                
                    new_gff.write(new_gff_string)
                    full_gff += new_gff_string
                exonless_contig = remove_exons("intronerate.gff","sequences/intron/{}_supercontig.fasta".format(gene))
                SeqIO.write(exonless_contig,"sequences/intron/{}_introns.fasta".format(gene),'fasta')
                
                os.chdir(basedir)
            else:
                sys.stderr.write("ERROR: Files not found for {} gene {}!\n".format(prefix,gene))
    with open("{}_genes.gff".format(prefix),'w') as all_gff:
        all_gff.write(full_gff)        



if __name__ == "__main__":main()
