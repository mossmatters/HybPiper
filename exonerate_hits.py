#!/usr/bin/env python

# Given an assembly file of genomic DNA and a file containing target proteins:
#     1. Use exonerate to determine hits for each contig.
#     2. Load contigs into a Biopython seqRecord dictionary.
#     3. For each protein hit, create a FASTA file with all contigs.
#         a. The contigs will be in order of their hit location to the target proteins.

######REQUIREMENTS########
# Python 2.6 or later
# exonerate in your $PATH
# Biopython
##########################

import sys, os, subprocess,math,argparse,logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#id_threshold = 55 #Percent identity between protein and contig hit.
first_search_filename = "exonerate_results.fasta"

def initial_exonerate(proteinfilename, assemblyfilename,prefix):
    """Conduct exonerate search, returns a dictionary of results.
    Using the ryo option in exonerate, the header should contain all the useful information."""
    logger = logging.getLogger("pipeline")
    
    outputfilename = "%s/exonerate_results.fasta" %prefix
    exonerate_ryo = '">%ti,%qi,%qab,%qae,%pi,(%tS),%tab,%tae\\n%tcs\\n"'
    exonerate_command = "exonerate -m protein2genome --showalignment no --showvulgar no -V 0 --ryo %s %s %s >%s" % (exonerate_ryo,proteinfilename,assemblyfilename,outputfilename)
    
    logger.debug(exonerate_command)
    #print exonerate_ryo
    #proc = subprocess.Popen(['exonerate','-m','protein2genome','--showalignment','no','-V','0','--showvulgar','no','--ryo',exonerate_ryo,proteinfilename,assemblyfilename])
    proc = subprocess.call(exonerate_command,shell=True)
    protHitsCount = 0
    #proc.wait()
    records = SeqIO.to_dict(SeqIO.parse(outputfilename,'fasta'))
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

def sort_key(elem):
    '''Sort by start location (increasing) then by end location (increasing), then by depth (decreasing)'''
    return elem[0],elem[1],-elem[2]

def get_contig_order(prot):
    """Given the dictionary of hits for a protein, return the dictionary with the fields sorted by start location."""
    logger = logging.getLogger("pipeline")
    
    tuplist =[(prot["hit_start"][i],prot["hit_end"][i],float(prot["assemblyHits"][i].split(",")[0].split("_")[5])) for i in range(len(prot["hit_start"]))]
    logger.debug("before sorting: {}".format(" ".join(prot["assemblyHits"])))
    logger.debug( tuplist )
    sorting_order = sorted(list(range(len(tuplist))),key=lambda k:sort_key(tuplist[k]))
    
    prot["assemblyHits"] = [prot["assemblyHits"][i] for i in sorting_order]
    prot["hit_start"] = [prot["hit_start"][i] for i in sorting_order]
    prot["hit_end"] = [prot["hit_end"][i] for i in sorting_order]
    prot["percentid"] = [prot["percentid"][i] for i in sorting_order]
    prot["hit_strand"] = [prot["hit_strand"][i] for i in sorting_order]
    
    logger.debug("After sorting: {}".format(" ".join(prot["assemblyHits"])))    
    return prot

def filter_by_percentid(prot,thresh):
    """Given a protein dictionary, return a protein dictionary minus entries with percentID below a threshold"""
    kept_indicies = [i for i in range(len(prot["percentid"])) if prot["percentid"][i] > thresh]
    return keep_indicies(kept_indicies,prot)        

def supercontig_exonerate(supercontig,protseq,prefix,thresh=55):
    """Given a long, joined contig and a protein sequence, return the exonerate hit(s)"""
    logger = logging.getLogger("pipeline")

    exonerate_ryo = '>%ti,%qi,%qab,%qae,%pi,(%tS)\\n%tcs\\n'
    temp_prot_filename = "%s/temp.prot.fa"%prefix
    temp_contig_filename =  "%s/temp.contig.fa"%prefix
    SeqIO.write(protseq,temp_prot_filename,'fasta')
    SeqIO.write(supercontig,temp_contig_filename,'fasta')
    logger.debug("Conducting exonerate search on supercontig")
    proc = subprocess.Popen(['exonerate','-m','protein2genome','--showalignment','no','-V','0','--showvulgar','no','--ryo',exonerate_ryo,temp_prot_filename,temp_contig_filename],stdout=subprocess.PIPE,universal_newlines=True)

    proc.wait()
    #print proc.stdout.readlines()
    supercontig_cds = [i for i in SeqIO.parse(proc.stdout,'fasta') if float(i.id.split(",")[4])>thresh]
    logger.debug("Supercontig lengths: %s" % " ".join([str(len(x.seq)) for x in supercontig_cds]))
    return supercontig_cds

def sort_byhitloc(seqrecord):
    """Key function for sorting based on the start location of a hit record."""
    return int(seqrecord.id.split(",")[2])

def subsume_supercontigs(supercontigs):
    """If one supercontig has a start and end location greater than all the others, throw the rest out"""
    logger = logging.getLogger("pipeline")
    supercontig_rangelist = [(int(x.id.split(",")[2]),int(x.id.split(",")[3])) for x in supercontigs]
    supercontig_ids = [x.id for x in supercontigs]
    logger.debug("Checking these ranges for supercontig: ")
    logger.debug(supercontig_rangelist)
    seqs_to_keep = range_connectivity(supercontig_rangelist,supercontig_ids)
    logger.debug("Keeping these contigs: ")
    logger.debug([supercontigs[x].id for x in seqs_to_keep])
    return [supercontigs[x] for x in seqs_to_keep]

def write_exonerate_stats(contig_id_list,prefix):
    '''Given a list of IDs from initial exonerate search, write info to a standard file'''
    with open("{}/exonerate_stats.csv".format(prefix),'w') as exonerate_statsfile:
        exonerate_statsfile.write("\n".join(contig_id_list)+'\n')

    
def fullContigs(prot,sequence_dict,assembly_dict,protein_dict,prefix,thresh=55):
    """Generates a contig from all hits to a protein. 
    If more than one hit, conduct a second exonerate search with the original contigs
    stitched together."""
    logger = logging.getLogger("pipeline")
    #logger.setLevel(logger.debug)
    numHits = len(prot["assemblyHits"])
    sequence_list = []
    contigHits = []
    
    logger.debug("All hits:")
    logger.debug(prot["assemblyHits"])
    write_exonerate_stats(prot["assemblyHits"],prefix)

    
    #print numHits
    if numHits == 1:
        #if prot["hit_strand"][0] == "+":
        return str(sequence_dict[prot["assemblyHits"][0]].seq)    #If only one hit to this protein.
        #else:
        #    return str(sequence_dict[prot["assemblyHits"][0]].seq.reverse_complement())
    else:
        for hit in range(len(prot["assemblyHits"])):
            assembly_seq_name = prot["assemblyHits"][hit].split(",")[0]
            logger.debug("Protein hit {} from {} to {} with {}% id on strand {}".format(assembly_seq_name,
                         prot["hit_start"][hit],
                         prot["hit_end"][hit],
                         prot["percentid"][hit],
                         prot["hit_strand"][hit]
                         ))
            if assembly_seq_name not in contigHits:        #Only add each contig once.
                if prot["hit_strand"][hit] == "+":
                    sequence_list.append(assembly_dict[assembly_seq_name])
                else:
                    sequence_list.append(SeqRecord(assembly_dict[assembly_seq_name].reverse_complement().seq,id=assembly_seq_name))
                contigHits.append(assembly_seq_name)
        #print("assembly_dict in fullContigs:")
        #print(assembly_dict)
        logger.debug("Contig order: {}".format(",".join([x.id for x in sequence_list])))
        logger.debug(",".join(contigHits))
#        logger.debug(assembly_dict.keys())
#     logger.debug([i for i in prot["assemblyHits"]])
#    logger.debug([(prot["hit_start"][i],prot["hit_end"][i]) for i in range(len(prot["hit_start"]))])
#     logger.debug(prot["hit_strand"])
#     logger.debug(prot["percentid"])
#     logger.debug("\n".join(["{}   {}".format(x,assembly_dict[x].seq) for x in contigHits]))
    supercontig = SeqRecord(Seq("".join(str(b.seq) for b in sequence_list)),id=prot["name"])
    logger.debug(">supercontig\n{}".format(supercontig.seq))
    #Need to remove contigs if they have the same basename
    supercontig_cds = supercontig_exonerate(supercontig,protein_dict[prot["name"]],prefix,thresh)
    if not supercontig_cds:
        sys.stderr.write("Supercontig below percent identity threshold!\n")
        return None
    logger.debug(" ".join(str(len(x)) for x in supercontig_cds))
    #Sort the supercontigs by hit location to the protein.
    joined_supercontig_cds = [b for b in supercontig_cds]
    joined_supercontig_cds.sort(key=sort_byhitloc)
    #print([x.id for x in joined_supercontig_cds])
    #logger.info([x for x in prot['assemblyHits'] if x in sequence_list])
    #write_exonerate_stats([x for x in prot['assemblyHits'] if x in sequence_list])

    #Get rid of supercontig sequences that are subsumed by longer sequences on the same stretch.
    joined_supercontig_cds = subsume_supercontigs(joined_supercontig_cds)
    
    
    SeqIO.write(joined_supercontig_cds,'%s/supercontig_exonerate.fasta'%prefix,'fasta') 
    if len(joined_supercontig_cds) == 1:
        logger.debug("One sequence remaining")
        return str(joined_supercontig_cds[0].seq)
    #One more Exonerate, just to be sure.
    superdupercontig = SeqRecord(Seq("".join(str(b.seq) for b in joined_supercontig_cds)),id=prot["name"])
    logger.debug(">joined_supercontig\n{}".format(superdupercontig.seq))
    #final_supercontig = [x for x in supercontig_exonerate(superdupercontig,protein_dict[prot["name"]],prefix)]
    #final_supercontig.sort(key=sort_byhitloc)
    #final_supercontig = subsume_supercontigs(final_supercontig)
    
    
    #return str(Seq("".join(str(b.seq) for b in final_supercontig)))
    return str(Seq("".join(str(b.seq) for b in joined_supercontig_cds)))        
    #print joined_supercontig_cds
    #print ""
    #return joined_supercontig_cds


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

def overlapping_contigs(prot,length_pct,depth_multiplier):
    """Given a protein dictionary, determine whether the hit ranges are overlapping,
    and save only those contigs that are not completely subsumed by other contigs."""
    logger = logging.getLogger("pipeline")
    range_list = [(prot["hit_start"][i],prot["hit_end"][i]) for i in range(len(prot["hit_start"]))]
    
    logger.debug(range_list)
    kept_indicies = range_connectivity(range_list,prot["assemblyHits"],prot_length = prot["reflength"],length_pct = length_pct,depth_multiplier=depth_multiplier)
    logger.debug(kept_indicies)
    return keep_indicies(kept_indicies,prot)


def best_by_percent_id(assemblyHits,full_length_indicies):
    '''Given a list of contig names, return the one with the best percent identity (fourth comma delimited field)'''
    logger = logging.getLogger("pipeline")
    max_percentid = 0
    for i in range(len(full_length_indicies)):
        percentid = float(assemblyHits[full_length_indicies[i]].split(",")[4])
        if percentid > max_percentid:
            logger.debug("percent_id: {}, maxpercent_id: {}".format(percentid,max_percentid))
            to_keep = full_length_indicies[i]
            max_percentid = percentid
    return to_keep
    
def best_by_depth(assemblyHits,full_length_indicies,thresh=10):
    '''If one contig has a depth that is 10x more than all the others, return that one, else return None'''
    logger=logging.getLogger("pipeline")
    depths = []
    for i in range(len(full_length_indicies)):
        depths.append((full_length_indicies[i],float(assemblyHits[full_length_indicies[i]].split(',')[0].split("_")[5])))
    depths.sort(reverse=True,key=lambda x: x[1])
    logger.debug(depths)
    depth_threshold = depths[0][1] / thresh
    logger.debug("Depth threshold: {}".format(depth_threshold))
    top_depth_best = all(i[1] <= depth_threshold for i in depths[1:]) 
    if top_depth_best:
        best_depth_contig = depths[0][0]
        logger.debug("Contig {} with depth {} is more than {} times greater depth than other contigs".format(best_depth_contig,depths[0][1],thresh))
        return best_depth_contig
    logger.debug("All contigs have similar depth")

    return None    


def range_connectivity(range_list,assemblyHits=None,prot_length=None,length_pct = 1,depth_multiplier = None,use_depth=False):
    """Given two sorted lists, representing the beginning and end of a range,
    Determine "connectivity" between consecutive elements of the list.
    For each connected segment, determine whether one segement "subsumes" the other."""
    
    logger = logging.getLogger("pipeline")

    starts = [a[0] for a in range_list]
    ends = [a[1] for a in range_list]
    
    if depth_multiplier:
        use_depth = True
    
    subsumed_ranges = []
    collapsed_ranges = []
    full_length_indicies = []
    num_breaks = 0
    if prot_length:
        max_length = prot_length
    else:
        max_length = max(ends) - min(starts)
    
    for i in range(len(range_list)):
            if abs(starts[i] - ends[i]) > max_length * length_pct:
                logger.debug("including long contig {}".format(range_list[i]))    
                full_length_indicies.append(i)
                subsumed_ranges = [range_list[i]]
            elif starts[i] == min(starts) and ends[i] == max(ends):
                logger.debug("Contig {} has range that subsumes all others!".format(i))
                subsumed_ranges = [range_list[i]]
                full_length_indicies.append(i)
            else:
                if len(full_length_indicies) > 0:
                    logger.debug("removing {}".format(range_list[i]))
                else:
                    subsumed_ranges.append(range_list[i])
    
    #If there are multiple full length hits, return the one with the best percent identity.
    if assemblyHits:
        if len(full_length_indicies) > 1:
            if use_depth:
                to_keep = best_by_depth(assemblyHits,full_length_indicies,depth_multiplier)
                if to_keep:
                    return [to_keep]
                else:
                    to_keep = best_by_percent_id(assemblyHits,full_length_indicies)        
                    return [to_keep]
            else:
                to_keep = best_by_percent_id(assemblyHits,full_length_indicies)        
                return [to_keep]
                                
    #If multiple contigs start at the same minimum (or end at the same maximum), keep the longest ones.
    subsumed_indices=[]
    if len(subsumed_ranges) > 1:
        logger.debug("SUBSUMING")
        best_start_end = 0
        best_end_start = 1000000000
        for i,r1 in enumerate(subsumed_ranges):
            for j,r2 in enumerate(subsumed_ranges):
                if i != j:
                    if tuple_subsume(r1,r2):
                        subsumed_indices.append(j)
        subsumed_set = set(subsumed_indices)
        kept_indices = [x for x in range(len(subsumed_ranges)) if x not in subsumed_set]
        return kept_indices
                    
#         for j in range(len(subsumed_ranges)):
#             if subsumed_ranges[j][0] == min(starts):
#                 if subsumed_ranges[j][1] > best_start_end:
#                     best_start_end = subsumed_ranges[j][1]
#                     longest_left = j
# 
#             elif subsumed_ranges[j][1] == max(ends):
#                 if subsumed_ranges[j][0] < best_end_start:
#                     best_end_start = subsumed_ranges[j][0]
#                     longest_right = j
#             else:
#                 collapsed_ranges.append(subsumed_ranges[j])
        
#         logger.debug("Best end start: {}".format(best_end_start))
#         logger.debug("Best start end: {}".format(best_start_end))
#         collapsed_ranges.append(subsumed_ranges[longest_left])
#         collapsed_ranges.append(subsumed_ranges[longest_right])
    else:
        collapsed_ranges = subsumed_ranges        

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
    
        
def tuple_overlap(a,b):
    """Given two tuples of length two, determine if the ranges overlap"""
    return a[0] < b[0] < a[1] or b[0] < a[0] < b[1]

def tuple_subsume(a,b):
    """Given two tuples of length two, determine if a has a range that includes b"""
    if b[0] >= a[0] and b[1] <= a[1]:
        return True
    else:
        return False

def reciprocal_best_hit(prot,proteinHits):
    """Given a protein dictionary and the dictionary of all protein dictionaries,
        Return the protein dictionary minus any contigs that have higher percentage hits to other proteins."""
    logger = logging.getLogger("pipeline")

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
                    logger.debug("%s %s" %(contig, full_contigname))
                    otherHit_idx = proteinHits[otherProt]["assemblyHits"].index(full_contigname)
                    
                    target_ranges = [sorted((prot["target_begin"][contig_idx],prot["target_end"][contig_idx])),sorted((proteinHits[otherProt]["target_begin"][otherHit_idx],proteinHits[otherProt]["target_end"][otherHit_idx]))]
                    logger.debug(repr(target_ranges))
                    #Check that the two contig hits have overlapping ranges.
                    if tuple_overlap(target_ranges[0],target_ranges[1]):                 
                        logger.debug("%s %s"%(repr(prot["percentid"][contig_idx]),repr(proteinHits[otherProt]["percentid"][otherHit_idx])))
                        if prot["percentid"][contig_idx] < proteinHits[otherProt]["percentid"][otherHit_idx]:
                            logger.debug("contig %s is a better hit to %s" %(contigname,otherProt))
                            maxProt = proteinHits[otherProt]["name"]
                    else:
                        logger.debug("ranges did not overlap")
        if maxProt == protname:
            kept_indicies.append(contig_idx)
    return keep_indicies(kept_indicies,prot)        

def paralog_test(exonerate_hits,prot,prefix):
    """Gives a warning if there are multiple hits of long length to the same protein"""
    logger = logging.getLogger("pipeline")
    protlength = len(prot)
    hitlengths = [abs(int(x.split(",")[2]) - int(x.split(",")[3])) for x in exonerate_hits["assemblyHits"]]
    logger.debug("protein length: {}".format(protlength))
    logger.debug("Hit lengths:")
    logger.debug(hitlengths)
    longhits = [x > 0.75*protlength for x in hitlengths]
    if sum(longhits) > 1:
        sys.stderr.write("WARNING: Multiple long-length exonerate hits for {}. Check for paralogs!\n".format(prot.id))
        with open("{}/paralog_warning.txt".format(prefix),'w') as pw:
            for hit in range(len(exonerate_hits["assemblyHits"])):
                if longhits[hit]:
                    pw.write(prot.id+ "\t"+exonerate_hits["assemblyHits"][hit] + "\n")    

def myTranslate(nucl):
    """Given a raw sequence of nucleotides, return raw sequence of amino acids."""
    #print nucl
    nucseq = Seq(nucl)
    #print nucseq
    aminoseq = nucseq.translate()
    return str(aminoseq)

def report_no_sequences(protname):
    sys.stderr.write("No valid sequences remain for {}!\n".format(protname))

def help():
    print("USAGE: python hybseq_pipeline.py proteinfile assemblyfile prefix")
    print("The program Exonerate must be in your $PATH.")
    print("You must have BioPython installed")
    print("A protein and a nucleotide directory will be created in the current directory with the prefix.")
    return    

def main(): 
    
    parser = argparse.ArgumentParser(description="exonerate_hits.py; Generate gene-by-gene protein and nucleotide files from Bait Capture Assembly")
    #parser.add_argument("-v", "--verbose",help="Report progress of pipeline to stdout",
    #    action="store_const",dest="loglevel",const=logging.INFO, default=logging.WARNING)
    parser.add_argument("--debug",help="Print debugging information for development testing.",
        action="store_true",dest="loglevel",default=False)
    parser.add_argument("proteinfile",help="FASTA file containing one 'bait' sequence per protein.")
    parser.add_argument("assemblyfile",help="FASTA file containing DNA sequence assembly.")
    parser.add_argument("--prefix",help="""Prefix for directory, files, and sequences generated from this assembly. 
            If not specified, will be extracted from assembly file name.""",default=None)
    parser.add_argument("--no_sequences",help="Do not generate protein and nucleotide sequence files.", action="store_true",default=False)
    parser.add_argument("--first_search_filename",help="Location of previously completed Exonerate results. Useful for testing.",default="no")
    parser.add_argument("-t","--threshold",help="Threshold for Percent Identity between contigs and proteins. default = 55%%",default=55,type=int)
    parser.add_argument("--length_pct",help="Include an exonerate hit if it is at least as long as X percentage of the reference protein length. Default = 100%%",default=90,type=int)
    parser.add_argument("--depth_multiplier",help="Accept any full-length hit if it has a coverage depth X times the next best hit. Set to zero to not use depth. Default = 10",default=10,type=int)
    
    
    args = parser.parse_args()

    proteinfilename = args.proteinfile
    assemblyfilename = args.assemblyfile
    if args.prefix:
        prefix = args.prefix
        if os.path.exists(prefix):
            pass
        else:
            os.mkdir(prefix)
    else:
        prefix = os.path.basename(assemblyfilename).split(".")[0]    

    logger = logging.getLogger("pipeline")
    ch = logging.StreamHandler()
    logger.addHandler(ch)
    if args.loglevel:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    #formatter = logging.Formatter('[%(levelname)s] %(message)s') #Commenting this out stops messages from appearing twice.
    #handler = logging.StreamHandler()
    #handler.setFormatter(formatter)
    #logger.addHandler(handler)
    try:
        proteinfile = open(proteinfilename)
    except IOError:
        print("The file %s could not be opened!" %proteinfilename)
        return()
        
    try:
        assemblyfile = open(assemblyfilename)
    except IOError:
        print("The file %s could not be opened!" % assemblyfilename)
        return()
    assembly_dict = SeqIO.to_dict(SeqIO.parse(assemblyfile,'fasta'))
    protein_dict = SeqIO.to_dict(SeqIO.parse(proteinfile,'fasta'))
    
    if os.path.exists(args.first_search_filename):     #Shortcut for Testing purposes
        logger.info("Reading initial exonerate results from file {}.".format(first_search_filename))
        sequence_dict = SeqIO.to_dict(SeqIO.parse(first_search_filename,'fasta'))
    else:
        #logger.info("Starting exonerate search, please wait.")
        sequence_dict = initial_exonerate(proteinfilename,assemblyfilename,prefix)
    proteinHits = protein_sort(sequence_dict)

    sys.stderr.write("There were {} exonerate hits for {}.\n".format(len(sequence_dict),proteinfilename))
    #print "There were %i unique proteins hit." % len(proteinHits)
    
    directory_name = "%s/sequences/FNA" % prefix
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)

    directory_name = "%s/sequences/FAA" % prefix
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    
    for prot in proteinHits:
        sys.stderr.write(prot)
        logger.debug(prot)
         #Put contigs in order along the protein.
         #logging.info("Searching for best hit to protein: %s" % proteinHits[prot]["name"])
#         logger.debug("Initial hits: %s" % " ".join(proteinHits[prot]["assemblyHits"]))
        logger.debug("Initial hits: %s" % len(proteinHits[prot]["assemblyHits"]))

        paralog_test(proteinHits[prot],protein_dict[prot],prefix)
        proteinHits[prot]["reflength"] = len(protein_dict[prot])
        
        proteinHits[prot] = get_contig_order(proteinHits[prot])
#         logger.debug("After get_contig_order: %s" % " ".join(proteinHits[prot]["assemblyHits"]))
        logger.debug("After get_contig_order: %d" % len(proteinHits[prot]["assemblyHits"]))
        #Remove contigs that are suboptimal hits. Only one protein hit allowed per contig.

        proteinHits[prot] = reciprocal_best_hit(proteinHits[prot],proteinHits)
#         logger.debug("After RBH: %s" % " ".join(proteinHits[prot]["assemblyHits"]))
        logger.debug("After RBH: %d" % len(proteinHits[prot]["assemblyHits"]))
        if len(proteinHits[prot]["assemblyHits"]) == 0:
            report_no_sequences(proteinHits[prot]["name"])
            continue        #All hits have been filtered out
         
         #Filter out contigs with a hit below a threshold
        proteinHits[prot] = filter_by_percentid(proteinHits[prot],args.threshold)
#         logger.debug("After filter_by_percent_id: %s" % " ".join(proteinHits[prot]["assemblyHits"]))
        logger.debug("After filter_by_percent_id: %d" % len(proteinHits[prot]["assemblyHits"]))
        if len(proteinHits[prot]["assemblyHits"]) == 0:
            report_no_sequences(proteinHits[prot]["name"])
            continue        #All hits have been filtered out
         
         #Delete contigs if their range is completely subsumed by another hit's range.
        proteinHits[prot] = overlapping_contigs(proteinHits[prot],args.length_pct*0.01,args.depth_multiplier)
#         logger.debug("After overlapping_contigs: %s" % " ".join(proteinHits[prot]["assemblyHits"]))
        logger.debug("After overlapping_contigs: %d" % len(proteinHits[prot]["assemblyHits"]))
         #Stitch together a "supercontig" containing all the hits and conduct a second exonerate search.    
        if len(proteinHits[prot]["assemblyHits"]) == 0:
            report_no_sequences(proteinHits[prot]["name"])
            continue        #All hits have been filtered out
        #print("sequence_dict")
        #print(sequence_dict)
        #print("assembly_dict:")
        #print(assembly_dict)
        nucl_sequence = fullContigs(proteinHits[prot],sequence_dict,assembly_dict,protein_dict,prefix,args.threshold)
        if nucl_sequence:
            if args.no_sequences:
                continue
            else:
                amino_sequence = myTranslate(nucl_sequence)
                seqID = prefix.split("/")[-1].strip("/")
                sys.stderr.write("Writing amino acid sequence, length: {}\n".format(len(amino_sequence)))
                sys.stdout.write("{}\t{}\n".format(prot.split("-")[-1],len(amino_sequence)))
                amino_filename = "%s/sequences/FAA/%s.FAA" % (prefix,prot.split("-")[-1])
                amino_file = open(amino_filename,'w')
                amino_file.write(">%s\n%s\n" % (seqID,amino_sequence))
                amino_file.close()
        
                nucleo_filename = "%s/sequences/FNA/%s.FNA" % (prefix,prot.split("-")[-1])
                nucleo_file = open(nucleo_filename,'w')
                nucleo_file.write(">%s\n%s\n" % (seqID,nucl_sequence))
                nucleo_file.close()
#      if "temp.contig.fa" in os.listdir(prefix):    
#         os.remove("%s/temp.contig.fa" % prefix)
#         os.remove("%s/temp.prot.fa" % prefix)
    proteinfile.close()
    assemblyfile.close()
    

if __name__ == "__main__":main()