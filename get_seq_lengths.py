
import os,sys
from Bio import SeqIO
#outfile = open("chloroplast_lengths.txt",'w')

helptext = '''

Usage: python baitfile.fasta namelist.txt dna/aa

    Prepare a list of names from each of your samples that used the HybSeqPipeline, one name per line.
    Indicate whether the bait file is DNA or AA (affects length calculations).
    
    The script prints a table to stdout. The first line contains gene names.
    The second line contains the length of the reference sequences (baits). 
    If there are multiple baits per gene, the mean length is reported.
    All other rows contain one sample per line.
    
    This script requires BioPython to parse the bait FASTA file.
    '''

if len(sys.argv) < 4:
    print(helptext)
    sys.exit()

baitfile = sys.argv[1]
namelistfile = sys.argv[2]
sequenceType = sys.argv[3]

if sequenceType.upper() == 'DNA':
    filetype = 'FNA'
elif sequenceType.upper() == 'AA':
    filetype = 'FAA'
else:
    print(helptext)
    sys.exit()

if not os.path.isfile(baitfile):
    print("Baitfile {} not found!".format(baitfile))
    sys.exit()
    
if not os.path.isfile(namelistfile):
    print("Name list file {} not found!".format(namelistfile))
    sys.exit()    

namelist =  [n.rstrip() for n in open(namelistfile).readlines()]

gene_names = []
reference_lengths = {}
for prot in SeqIO.parse(baitfile,"fasta"):
    protname = prot.id.split("-")[-1]
    gene_names.append(protname)
    if protname in reference_lengths:
        reference_lengths[protname].append(len(prot.seq))
    else:
        reference_lengths[protname] = [len(prot.seq)]
unique_names = list(set(gene_names))
avg_ref_lengths = [(sum(reference_lengths[gene])/len(reference_lengths[gene])) for gene in unique_names]
sys.stdout.write("Species\t{}\nMeanLength\t{}\n".format("\t".join(unique_names),"\t".join([str(x) for x in avg_ref_lengths])))

for name in namelist:
    parentDir,name = os.path.split(name)
    if not name:
        parentDir,name = os.path.split(parentDir)
    name_lengths = []
    for gene in range(len(unique_names)):
        read_file = os.path.join(parentDir,name,unique_names[gene],name,"sequences",filetype,"{}.{}".format(unique_names[gene],filetype))
        if os.path.exists(read_file):
            seq_length = len(SeqIO.read(read_file,'fasta').seq)
            if seq_length > 1.5 * avg_ref_lengths[gene]:
                sys.stderr.write("****WARNING! Sequence length for {} is more than 50% longer than {} reference!\n".format(name,unique_names[gene]))
            name_lengths.append(str(seq_length))
        else:
            name_lengths.append("0")
        
    sys.stdout.write("{}\t{}\n".format(name,"\t".join(name_lengths)))

