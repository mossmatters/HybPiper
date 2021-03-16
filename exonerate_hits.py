#!/usr/bin/env python

################################################################################################
#### CJJ NOTE THAT THIS SCRIPT HAS BEEN HEAVILY ALTERED TO ADD THE FOLLOWING FUNCTIONALITY: ####
################################################################################################

"""
- When stitching together contigs tos create a supercontig, trim the Exonerate hits to hit sequence only, rather than
  using the whole contig sequence,.
- For supercontigs, trim overlaps between Exonerate hits when stitching together sequences. To me, these overlaps seem
  most likely to occur when when Exonerate hits come from partially assembled paralogs, and the contigs have
  overlapping termini (which presumably end when the SPAdes assembly graphs get too confused... ).
"""



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

import sys, os, subprocess, math, argparse, logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import tee  # CJJ
import re  # CJJ

# id_threshold = 55 #Percent identity between protein and contig hit.
first_search_filename = "exonerate_results.fasta"


def file_exists_and_not_empty(file_name):  # CJJ added for initial exonerate tests
    """
    Check if file exists and is not empty by confirming that its size is not 0 bytes
    """
    # Check if file exist and is not empty
    return os.path.isfile(file_name) and not os.path.getsize(file_name) == 0


def initial_exonerate(proteinfilename, assemblyfilename, prefix):
    """Conduct exonerate search, returns a dictionary of results.
    Using the ryo option in exonerate, the header should contain all the useful information."""
    logger = logging.getLogger("pipeline")
    outputfilename = "%s/exonerate_results.fasta" % prefix
    exonerate_ryo = '">%ti,%qi,%qab,%qae,%pi,(%tS),%tab,%tae\\n%tcs\\n"'

    exonerate_command = "exonerate -m protein2genome --showalignment no --showvulgar no -V 0 --refine full --ryo %s " \
                        "%s %s >%s" % (exonerate_ryo, proteinfilename, assemblyfilename, outputfilename)  # CJJ refine
    logger.debug(exonerate_command)
    try:
        proc = subprocess.call(exonerate_command, shell=True)
    except:
        logger.debug(f'exonerate_command with refine failed for {outputfilename}')
    if file_exists_and_not_empty(outputfilename):  # CJJ
        records = SeqIO.to_dict(SeqIO.parse(outputfilename, 'fasta'))
        return records
    else:
        exonerate_command = "exonerate -m protein2genome --showalignment no --showvulgar no -V 0 --ryo %s %s %s >%s" \
                            % (exonerate_ryo, proteinfilename, assemblyfilename, outputfilename)
        logger.debug(exonerate_command)
        proc = subprocess.call(exonerate_command, shell=True)
        records = SeqIO.to_dict(SeqIO.parse(outputfilename, 'fasta'))
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
            proteinHits[protein] = {"assemblyHits": [",".join(hit)],
                                    "hit_start": [int(hit[2])],
                                    "hit_end": [int(hit[3])],
                                    "percentid": [float(hit[4])],
                                    "hit_strand": [hit[5][1]],
                                    "target_begin": [int(hit[6])],
                                    "target_end": [int(hit[7])],
                                    "name": protein
                                    }
    return proteinHits


def sort_key(elem):
    """Sort by start location (increasing) then by end location (increasing), then by depth (decreasing)"""
    return elem[0], elem[1], -elem[2]


def get_contig_order(prot):
    """Given the dictionary of hits for a protein, return the dictionary with the fields sorted by start location."""
    logger = logging.getLogger("pipeline")
    tuplist = [(prot["hit_start"][i], prot["hit_end"][i], float(prot["assemblyHits"][i].split(",")[0].split("_")[5]))
               for i in range(len(prot["hit_start"]))]
    logger.debug("before sorting: {}".format(" ".join(prot["assemblyHits"])))
    logger.debug(tuplist)
    sorting_order = sorted(list(range(len(tuplist))), key=lambda k: sort_key(tuplist[k]))
    prot["assemblyHits"] = [prot["assemblyHits"][i] for i in sorting_order]
    prot["hit_start"] = [prot["hit_start"][i] for i in sorting_order]
    prot["hit_end"] = [prot["hit_end"][i] for i in sorting_order]
    prot["percentid"] = [prot["percentid"][i] for i in sorting_order]
    prot["hit_strand"] = [prot["hit_strand"][i] for i in sorting_order]
    prot["target_begin"] = [prot["target_begin"][i] for i in sorting_order]  # CJJ ADDED
    prot["target_end"] = [prot["target_end"][i] for i in sorting_order]  # CJJ ADDED
    logger.debug("After sorting: {}".format(" ".join(prot["assemblyHits"])))
    return prot


def filter_by_percentid(prot, thresh):
    """Given a protein dictionary, return a protein dictionary minus entries with percentID below a threshold"""
    kept_indicies = [i for i in range(len(prot["percentid"])) if prot["percentid"][i] > thresh]
    return keep_indicies(kept_indicies, prot)


def supercontig_exonerate(supercontig, protseq, prefix, thresh=55):
    """Given a long, joined contig and a protein sequence, return the exonerate hit(s)"""
    logger = logging.getLogger("pipeline")
    exonerate_ryo = '>%ti,%qi,%qab,%qae,%pi,(%tS)\\n%tcs\\n'
    temp_prot_filename = "%s/temp.prot.fa" % prefix
    temp_contig_filename = "%s/temp.contig.fa" % prefix
    SeqIO.write(protseq, temp_prot_filename, 'fasta')
    SeqIO.write(supercontig, temp_contig_filename, 'fasta')
    logger.debug("Conducting exonerate search on supercontig")
    proc = subprocess.Popen(
        ['exonerate', '-m', 'protein2genome', '--showalignment', 'no', '-V', '0', '--showvulgar', 'no', '--ryo',
         exonerate_ryo, temp_prot_filename, temp_contig_filename], stdout=subprocess.PIPE, universal_newlines=True)
    proc.wait()
    supercontig_cds = [i for i in SeqIO.parse(proc.stdout, 'fasta') if float(i.id.split(",")[4]) > thresh]
    logger.debug("Supercontig lengths: %s" % " ".join([str(len(x.seq)) for x in supercontig_cds]))
    return supercontig_cds


def sort_byhitloc(seqrecord):
    """Key function for sorting based on the start location of a hit record."""
    return int(seqrecord.id.split(",")[2])


def subsume_supercontigs(supercontigs):
    """If one supercontig has a start and end location greater than all the others, throw the rest out"""
    logger = logging.getLogger("pipeline")
    supercontig_rangelist = [(int(x.id.split(",")[2]), int(x.id.split(",")[3])) for x in supercontigs]
    supercontig_ids = [x.id for x in supercontigs]
    logger.debug("Checking these ranges for supercontig: ")
    logger.debug(supercontig_rangelist)
    seqs_to_keep = range_connectivity(supercontig_rangelist, supercontig_ids)
    logger.debug("Keeping these contigs: ")
    logger.debug([supercontigs[x].id for x in seqs_to_keep])
    return [supercontigs[x] for x in seqs_to_keep]


def write_exonerate_stats(contig_id_list, prefix):
    """Given a list of IDs from initial exonerate search, write info to a standard file"""
    with open("{}/exonerate_stats.csv".format(prefix), 'w') as exonerate_statsfile:
        exonerate_statsfile.write("\n".join(contig_id_list) + '\n')


def write_genes_with_supercontigs(data, prefix):  # CJJ
    """Write a file listing genes for which a supercontig was created. These per sample files are collated
    in the reads_first.py script after all genes have completed)."""
    with open("{}/genes_with_supercontigs.csv".format(prefix), 'w') as supercontig_reportfile:
        supercontig_reportfile.write(f'{data}\n')


def write_supercontigs_with_discordant_readpairs(data, prefix):  # CJJ
    """Write a file listing supercontigs for which one read maps perfectly and the other has mismatches with the
    reference. These per sample files are collated in reads_first.py script after all genes have
    completed)."""
    with open("{}/supercontigs_with_discordant_readpairs.csv".format(prefix), 'w') as discordant_supercontig_reportfile:
        discordant_supercontig_reportfile.write(f'{data}\n')


def pairwise(iterable):  # CJJ
    """s -> (s0,s1), (s1,s2), (s2, s3), ...
    Used in the function fullContigs to iterate over overlapping pairs of hit_start_and_end_indices.
    """
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def grouped(iterable, n):  # CJJ
    """s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ...
    Used in the function fullContigs to iterate over non-overlapping pairs of reads from a sam file (i.e. reads 1+2,
    then reads 3+4 etc).
    """
    return zip(*[iter(iterable)]*n)


def fullContigs(prot, sequence_dict, assembly_dict, protein_dict, prefix, thresh=55, nosupercontigs=False,
                interleaved_reads='None', memory=1, discordant_cutoff=100, edit_distance=7, threads=1):
    """Generates a contig from all hits to a protein.
    If more than one hit, conduct a second exonerate search with the original contigs
    stitched together."""

    sys.stderr.write(f'\nCJJ from within fullContigs function\n')
    sys.stderr.flush()

    logger = logging.getLogger("pipeline")
    numHits = len(prot["assemblyHits"])
    sequence_list = []
    contigHits = []

    logger.debug("All hits:")
    logger.debug(prot["assemblyHits"])
    write_exonerate_stats(prot["assemblyHits"], prefix)
    log_sample_and_gene = f'{os.path.basename(os.getcwd())},{prot["name"].split("-")[1]}'  # CJJ

    if numHits == 1:
        log_entry = f'{log_sample_and_gene},No supercontig created'  # CJJ
        write_genes_with_supercontigs(log_entry, prefix)
        return str(sequence_dict[prot["assemblyHits"][0]].seq)  # If only one hit to this protein.
    elif nosupercontigs:  # CJJ Recover the sequence from the longest single hit only.
        longest_hit_length = 0
        longest_hit_name = 'NA'
        hit_index = 0
        for hit in range(len(prot["assemblyHits"])):
            assembly_seq_name = prot["assemblyHits"][hit].split(",")[0]
            start = prot["hit_start"][hit]
            end = prot["hit_end"][hit]

            length = end - start
            if length >= longest_hit_length:
                longest_hit_length = length
                longest_hit_name = assembly_seq_name
                hit_index = hit

        log_entry = [log_sample_and_gene, 'Supercontig step skipped']
        write_genes_with_supercontigs(','.join(log_entry), prefix)
        return str(sequence_dict[prot["assemblyHits"][hit_index]].seq)
    else:
       ################################## CJJ trim supercontigs ########################################################

        # CJJ this is very verbose and convoluted, partly because I have to fit in with how HybPiper was written. But
        # there's no doubt a more elegant way of doing this. Review...eventually...
        hit_start_and_end_indices = []
        nucleotide_slice_indices = {}
        for hit in range(len(prot["assemblyHits"])):
            assembly_seq_name = prot["assemblyHits"][hit].split(",")[0]
            assembly_seq_length = len(assembly_dict[assembly_seq_name].seq)
            hit_start_and_end_indices.append([assembly_seq_name,
                                              prot["hit_start"][hit],
                                              prot["hit_end"][hit],
                                              prot["target_begin"][hit],
                                              prot["target_end"][hit],
                                              prot["hit_strand"][hit],
                                              assembly_seq_length])
            nucleotide_slice_indices[assembly_seq_name] = [None, None]  # CJJ set default slice indices
        # print(hit_start_and_end_indices)
        for pairs in pairwise(hit_start_and_end_indices):
            left_contig_name = pairs[0][0]
            right_contig_name = pairs[1][0]
            left_contig_strand = pairs[0][5]
            right_contig_strand = pairs[1][5]
            left_contig_length = pairs[0][6]
            right_contig_length = pairs[1][6]


            left_prot_target_start = int(pairs[0][1])
            left_prot_target_end = int(pairs[0][2])
            # print(f'left_prot_target_end: {left_prot_target_end}')
            right_prot_target_start = int(pairs[1][1])
            # print(f'right_prot_target_start: {right_prot_target_start}')

            if left_prot_target_end > right_prot_target_start:  # Check for overlap between the Exonerate hits
                overlap_offset_in_nucleotides = (left_prot_target_end - right_prot_target_start) * 3
                # print(f'overlap_offset_in_nucleotides: {overlap_offset_in_nucleotides}')

                ########################################################################################################
                ########################################################################################################
                # Calculate intron offset:

                intron_offset = 0 # set as default for slice calculations below
                left_query_span_nucleotides_in_exonerate_hit = (left_prot_target_end - left_prot_target_start) * 3 # i.e. how many
                # nucleotides in the contig query are actually in the Exonerate hit (so, not counting introns).
                # print(f'left_query_span_nucleotides_in_exonerate_hit: {left_query_span_nucleotides_in_exonerate_hit}')

                if left_contig_strand == '+':
                    left_query_spades_contig_start = int(pairs[0][3])  # i.e. smaller number relative to spades contig
                    left_query_spades_contig_end = int(pairs[0][4])  # i.e. larger number relative to spades contig
                    left_query_span_nucleotides_in_spades_contig = left_query_spades_contig_end - \
                                                                   left_query_spades_contig_start
                    # print(f'Positive strand left_query_span_nucleotides_in_spades_contig: '
                    #       f'{left_query_span_nucleotides_in_spades_contig}')
                elif left_contig_strand == '-':
                    left_query_spades_contig_start = int(pairs[0][4]) # i.e. smaller number relative to spades contig
                    left_query_spades_contig_end = int(pairs[0][3])  # i.e. larger number relative to spades contig
                    left_query_span_nucleotides_in_spades_contig = left_query_spades_contig_end - \
                                                                   left_query_spades_contig_start
                    # print(f'Negative strand left_query_span_nucleotides_in_spades_contig: '
                    #       f'{left_query_span_nucleotides_in_spades_contig}')

                # Check if the left Exonerate hit contains introns:  CHANGE THIS TO CHECK IF THE INTRON IS IN THE
                # OVERLAP, OTHERWISE THE INTRON OFFSET SHOULDN'T BE USED!!! CAN I EVEN DO THIS WITH THE EXONERATE
                # RESULTS HYBPIPER GENERATES?
                if left_query_span_nucleotides_in_exonerate_hit < left_query_span_nucleotides_in_spades_contig:
                    # print('Query span does not equal target span!')
                    intron_offset = left_query_span_nucleotides_in_spades_contig - left_query_span_nucleotides_in_exonerate_hit
                    # print(f'Intron offset is: {intron_offset}')
                ########################################################################################################
                ########################################################################################################

                # Calculate slices:

                # Left hit:
                if left_contig_strand == '+': # CJJ get slice indexes if the left hit is on the positive strand
                    # print('yep left positive')
                    left_slice_start = int(pairs[0][3])  # i.e. the nucleotide associated with the left side of exonerate results, NOT from the contig
                    left_contig_hit_end = int(pairs[0][4])  # i.e. the nucleotide associated with the right side of exonerate results, NOT from the contig
                    # left_slice_end = int(pairs[0][4]) - overlap_offset_in_nucleotides - intron_offset
                    left_slice_end = int(left_contig_hit_end - overlap_offset_in_nucleotides)  # CJJ 13Oct2020

                elif left_contig_strand == '-':  # CJJ get slice indexes if the left hit is on the negative strand
                    left_contig_hit_start = int(pairs[0][3])  # Note that here hit_start mean relative to the translated protein target.
                    left_contig_hit_end = int(pairs[0][4])
                    left_slice_start = left_contig_length - left_contig_hit_start
                    # print(f'left_contig_strand negative; left_slice_start: {left_slice_start}')
                    # left_slice_end = left_slice_start + left_contig_hit_start - left_contig_hit_end - \
                    #                  overlap_offset_in_nucleotides - intron_offset
                    left_slice_end = left_slice_start + left_contig_hit_start - left_contig_hit_end - \
                                     overlap_offset_in_nucleotides
                    # print(f'left_contig_strand negative; left_slice_end: {left_slice_end}')

                # Right hit:
                if right_contig_strand == '+':
                    right_slice_start = int(pairs[1][3])
                elif right_contig_strand == '-':
                    right_contig_hit_start = int(pairs[1][3])  # Note that here hit_start mean relative to the translated protein target.
                    # print(f'right_contig_hit_start: {right_contig_hit_start}')
                    right_slice_start = right_contig_length - right_contig_hit_start  # i.e. trim off any contig sequence that isn't in Exonerate hit.
                    # print(f'right_slice_start; {right_slice_start}')


                ########################################################################################################
                ########################################################################################################

                # Populate slice dictionary for contig pair:

                if nucleotide_slice_indices[left_contig_name][0] is not None:
                    pass
                else:
                    nucleotide_slice_indices[left_contig_name][0] = left_slice_start
                if nucleotide_slice_indices[left_contig_name][1] is not None:
                    pass
                else:
                    nucleotide_slice_indices[left_contig_name][1] = left_slice_end
                if nucleotide_slice_indices[right_contig_name][0] is not None:
                    pass
                else:
                    nucleotide_slice_indices[right_contig_name][0] = right_slice_start

        # print(nucleotide_slice_indices)

        # Check if trimming has been performed and, if so, write slice indices to the log file
        trimming_performed = False
        for key, value in nucleotide_slice_indices.items():
            if value[0] is not None or value[1] is not None:
                trimming_performed = True
                break
        if trimming_performed:
            log_entry = [log_sample_and_gene]
            for key, value in nucleotide_slice_indices.items():
                value = re.sub(',', '', str(value))
                log_entry.append(f'{key} {value}')
        else:
            log_entry = [log_sample_and_gene, 'Supercontig created but no contig trimming performed']
        write_genes_with_supercontigs(','.join(log_entry), prefix)

        ############################ CJJ Finished getting slice indexes for trimming ###################################

        for hit in range(len(prot["assemblyHits"])):
            assembly_seq_name = prot["assemblyHits"][hit].split(",")[0]
            logger.debug("Protein hit {} from {} to {} with {}% id on strand {}".format(assembly_seq_name,
                                                                                        prot["hit_start"][hit],
                                                                                        prot["hit_end"][hit],
                                                                                        prot["percentid"][hit],
                                                                                        prot["hit_strand"][hit]
                                                                                        ))
            if assembly_seq_name not in contigHits:  # Only add each contig once. # CJJ: contigHits is an empty list
                # defined at the start of the function
                start = nucleotide_slice_indices[assembly_seq_name][0]
                end = nucleotide_slice_indices[assembly_seq_name][1]
                # print(f'start, end: {start}, {end}')
                if prot["hit_strand"][hit] == "+":
                    # sequence_list is an empty list defined at the start of the function
                    sequence_list.append(assembly_dict[assembly_seq_name][start:end])  # CJJ slice seq using indices
                else:
                    sequence_list.append(
                        SeqRecord(assembly_dict[assembly_seq_name].reverse_complement().seq[start:end],
                                  id=assembly_seq_name))
                contigHits.append(assembly_seq_name) # contigHits is an empty list defined at the start of the function
        logger.debug("Contig order: {}".format(",".join([x.id for x in sequence_list])))
        logger.debug(",".join(contigHits))
    supercontig = SeqRecord(Seq("".join(str(b.seq) for b in sequence_list)), id=prot["name"])
    logger.debug(">supercontig\n{}".format(supercontig.seq))
    # Need to remove contigs if they have the same basename

    ####################################################################################################################
    # Run Exonerate for a second time using the same protein query and the supercontig as subject
    ####################################################################################################################
    supercontig_cds = supercontig_exonerate(supercontig, protein_dict[prot["name"]], prefix, thresh)  # Returns a
    # list. Can return multiple fasta hits.
    if not supercontig_cds:
        sys.stderr.write("Supercontig below percent identity threshold!\n")
        return None
    logger.debug(" ".join(str(len(x)) for x in supercontig_cds))
    # Sort the supercontigs by hit location to the protein.
    joined_supercontig_cds = [b for b in supercontig_cds]  # CJJ: does this do anything? 07March2021
    joined_supercontig_cds.sort(key=sort_byhitloc)
    # Get rid of supercontig sequences that are subsumed by longer sequences on the same stretch.
    joined_supercontig_cds = subsume_supercontigs(joined_supercontig_cds)
    SeqIO.write(joined_supercontig_cds, '%s/supercontig_exonerate.fasta' % prefix, 'fasta') # This can be multiple
    # fasta seqs written to the same file.
    discordant_reads = 0
    discordant_cutoff = discordant_cutoff    # CJJ user modifiable
    edit_distance = edit_distance            # CJJ user modifiable
    maxindel = 0
    minid = 0.76                             # CJJ should be user modifiable
    if len(joined_supercontig_cds) == 1:  # CJJ: i.e. if file supercontig_exonerate.fasta contains a single fast seq.

        ################ CJJ mapping check 1: if only one sequence left after filtering above ##########################
        supercontig_reference = joined_supercontig_cds
        SeqIO.write(supercontig_reference, '%s/CJJ_supercontig.fasta' % prefix, 'fasta')
        logger.debug("One sequence remaining")
        sys.stderr.write(f'\nInterleaved_reads: {interleaved_reads}\n')
        sys.stderr.write(f'\nSupercontig_reference: {supercontig_reference}\n')
        sys.stderr.flush()

        # CJJ How to specify threads? This script is launched via parallel, so if I specify multiple threads here then
        # it'll overburden the number of cpus requested by the slurm job.
        bbmap_command = f'bbmap.sh -Xmx{memory}g -t={threads} ref={prefix}/CJJ_supercontig.fasta in={interleaved_reads} ' \
                        f'out={prefix}/CJJ_supercontig.sam interleaved=t pairedonly=t mappedonly=t ' \
                        f'maxindel={maxindel} strictmaxindel=t nodisk=t minid={minid} ambiguous=toss 2> /dev/null'
        # sys.stderr.write(f'\nbbmap_command: {bbmap_command}\n')
        # sys.stderr.flush()
        exitcode = subprocess.call(bbmap_command, shell=True)

        samtools_proper_pair = f'samtools view -f 3 -h {prefix}/CJJ_supercontig.sam > ' \
                               f'{prefix}/CJJ_supercontig_properPair.sam'
        # sys.stderr.write(f'\nbwa samtools_proper_pair command: {samtools_proper_pair}\n')
        # sys.stderr.flush()
        exitcode = subprocess.call(samtools_proper_pair, shell=True)
        samfile_reads = []
        with open(f'{prefix}/CJJ_supercontig_properPair.sam') as samfile:
            lines = samfile.readlines()
            for line in lines:
                if not line.startswith('@'):
                    samfile_reads.append(line)

        with open(f'{prefix}/diagnostic_reads.sam', 'w') as diagnostic_reads:
            for forward, reverse in grouped(samfile_reads, 2):
                forward_edit_distance = (forward.split('\t')[11]).split(':')[2]
                reverse_edit_distance = (reverse.split('\t')[11]).split(':')[2]
                if int(forward_edit_distance) == 0 and int(reverse_edit_distance) >= edit_distance:
                    discordant_reads += 1
                    diagnostic_reads.write(forward)
                    diagnostic_reads.write(reverse)
                elif int(reverse_edit_distance) == 0 and int(forward_edit_distance) >= edit_distance:
                    discordant_reads += 1
                    diagnostic_reads.write(forward)
                    diagnostic_reads.write(reverse)

        if discordant_reads > discordant_cutoff:
            log_entry = f'{log_sample_and_gene},POTENTIAL MULTIPLE PARALOGS PRESENT IN SUPERCONTIG'
            write_supercontigs_with_discordant_readpairs(log_entry, prefix)

        return str(joined_supercontig_cds[0].seq)
    # One more Exonerate, just to be sure.
    superdupercontig = SeqRecord(Seq("".join(str(b.seq) for b in joined_supercontig_cds)), id=prot["name"])
    logger.debug(">joined_supercontig\n{}".format(superdupercontig.seq))


    ##############  CJJ mapping check 2: if multiple sequences left after filtering above  #############################
    supercontig_reference = superdupercontig
    SeqIO.write(supercontig_reference, '%s/CJJ_supercontig.fasta' % prefix, 'fasta')
    sys.stderr.write(f'\nInterleaved_reads: {interleaved_reads}\n')
    sys.stderr.write(f'\nSuperdupercontig_reference: {supercontig_reference}\n')
    sys.stderr.flush()

    bbmap_command = f'bbmap.sh -Xmx{memory}g  -t={threads} ref={prefix}/CJJ_supercontig.fasta in={interleaved_reads} ' \
                    f'out={prefix}/CJJ_supercontig.sam interleaved=t pairedonly=t mappedonly=t maxindel={maxindel} ' \
                    f'strictmaxindel=t nodisk=t minid={minid} ambiguous=toss 2> /dev/null'

    # sys.stderr.write(f'\nbbmap_command: {bbmap_command}\n')
    # sys.stderr.flush()
    exitcode = subprocess.call(bbmap_command, shell=True)

    samtools_proper_pair = f'samtools view -f 3 -h {prefix}/CJJ_supercontig.sam > ' \
                           f'{prefix}/CJJ_supercontig_properPair.sam'
    # sys.stderr.write(f'\nbwa samtools_proper_pair command: {samtools_proper_pair}\n')
    # sys.stderr.flush()
    exitcode = subprocess.call(samtools_proper_pair, shell=True)
    samfile_reads = []
    with open(f'{prefix}/CJJ_supercontig_properPair.sam') as samfile:
        lines = samfile.readlines()
        for line in lines:
            if not line.startswith('@'):
                samfile_reads.append(line)

    with open(f'{prefix}/diagnostic_reads.sam', 'w') as diagnostic_reads:
        for forward, reverse in grouped(samfile_reads, 2):
            forward_edit_distance = (forward.split('\t')[11]).split(':')[2]
            reverse_edit_distance = (reverse.split('\t')[11]).split(':')[2]
            if int(forward_edit_distance) == 0 and int(reverse_edit_distance) >= edit_distance:
                discordant_reads += 1
                diagnostic_reads.write(forward)
                diagnostic_reads.write(reverse)
            elif int(reverse_edit_distance) == 0 and int(forward_edit_distance) >= edit_distance:
                discordant_reads += 1
                diagnostic_reads.write(forward)
                diagnostic_reads.write(reverse)

    if discordant_reads > discordant_cutoff:
        log_entry = f'{log_sample_and_gene},POTENTIAL MULTIPLE PARALOGS PRESENT IN SUPERCONTIG'
        write_supercontigs_with_discordant_readpairs(log_entry, prefix)

    return str(Seq("".join(str(b.seq) for b in joined_supercontig_cds)))


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
    """Given a list of indicies to keep and a protein dictionary, return the dictionary with only the specified entries
    remaining"""

    assHit = []
    hitstart = []
    hitend = []
    percentid = []
    strands = []
    targetbegin = []
    targetend = []

    for a in kept_indicies:
        # print('CJJ -keep_indicies- a:', a)
        assHit.append(prot["assemblyHits"][a])
        # print('CJJ -keep_indicies- assHit:', assHit)
        hitstart.append(prot["hit_start"][a])
        # print('CJJ -keep_indicies- hitstart:', hitstart)
        hitend.append(prot["hit_end"][a])
        # print('CJJ -keep_indicies- hitend:', hitend)
        percentid.append(prot["percentid"][a])
        strands.append(prot["hit_strand"][a])
        targetbegin.append(prot["target_begin"][a])
        # print('CJJ -keep_indicies- targetbegin:', targetbegin)
        targetend.append(prot["target_end"][a])
        # print('CJJ -keep_indicies- targetend:', targetend)

    prot["assemblyHits"] = assHit
    prot["hit_start"] = hitstart
    prot["hit_end"] = hitend
    prot["percentid"] = percentid
    prot["hit_strand"] = strands
    prot["target_begin"] = targetbegin
    prot["target_end"] = targetend

    return prot


def overlapping_contigs(prot, length_pct, depth_multiplier):
    """Given a protein dictionary, determine whether the hit ranges are overlapping,
    and save only those contigs that are not completely subsumed by other contigs."""
    logger = logging.getLogger("pipeline")
    range_list = [(prot["hit_start"][i], prot["hit_end"][i]) for i in range(len(prot["hit_start"]))]

    logger.debug(range_list)
    kept_indicies = range_connectivity(range_list, prot["assemblyHits"], prot_length=prot["reflength"],
                                       length_pct=length_pct, depth_multiplier=depth_multiplier)
    logger.debug(kept_indicies)
    return keep_indicies(kept_indicies, prot)


def best_by_percent_id(assemblyHits, full_length_indicies):
    """Given a list of contig names, return the one with the best percent identity (fourth comma delimited field)"""
    logger = logging.getLogger("pipeline")
    max_percentid = 0
    for i in range(len(full_length_indicies)):
        percentid = float(assemblyHits[full_length_indicies[i]].split(",")[4])
        if percentid > max_percentid:
            logger.debug("percent_id: {}, maxpercent_id: {}".format(percentid, max_percentid))
            to_keep = full_length_indicies[i]
            max_percentid = percentid
    return to_keep


def best_by_depth(assemblyHits, full_length_indicies, thresh=10):
    """If one contig has a depth that is 10x more than all the others, return that one, else return None"""
    logger = logging.getLogger("pipeline")
    depths = []
    for i in range(len(full_length_indicies)):
        depths.append(
            (full_length_indicies[i], float(assemblyHits[full_length_indicies[i]].split(',')[0].split("_")[5])))
    depths.sort(reverse=True, key=lambda x: x[1])
    logger.debug(depths)
    depth_threshold = depths[0][1] / thresh
    logger.debug("Depth threshold: {}".format(depth_threshold))
    top_depth_best = all(i[1] <= depth_threshold for i in depths[1:])
    if top_depth_best:
        best_depth_contig = depths[0][0]
        logger.debug(
            "Contig {} with depth {} is more than {} times greater depth than other contigs".format(
                best_depth_contig, depths[0][1], thresh))
        return best_depth_contig
    logger.debug("All contigs have similar depth")
    return None


def range_connectivity(range_list, assemblyHits=None, prot_length=None, length_pct=1, depth_multiplier=None,
                       use_depth=False):
    """Given two sorted lists, representing the beginning and end of a range,
    Determine "connectivity" between consecutive elements of the list.
    For each connected segment, determine whether one segment "subsumes" the other."""

    logger = logging.getLogger("pipeline")

    starts = [a[0] for a in range_list]  # range_list e.g.: [(25, 108), (32, 66), (100, 201)]
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
            subsumed_ranges = [range_list[i]]  # CJJ Note this replaces rather than appends; this means it'll keep the
            # contig with this range, I think.
        elif starts[i] == min(starts) and ends[i] == max(ends):
            logger.debug("Contig {} has range that subsumes all others!".format(i))
            subsumed_ranges = [range_list[i]]
            full_length_indicies.append(i)
        else:
            if len(full_length_indicies) > 0:
                logger.debug("removing {}".format(range_list[i]))
            else:
                subsumed_ranges.append(range_list[i])  # CJJ Note this appends rather than replaces

    # If there are multiple full length hits, return the one with the best percent identity.
    if assemblyHits:
        if len(full_length_indicies) > 1:
            if use_depth:
                to_keep = best_by_depth(assemblyHits, full_length_indicies, depth_multiplier)
                if to_keep:
                    return [to_keep]
                else:
                    to_keep = best_by_percent_id(assemblyHits, full_length_indicies)
                    return [to_keep]
            else:
                to_keep = best_by_percent_id(assemblyHits, full_length_indicies)
                return [to_keep]

    # If multiple contigs start at the same minimum (or end at the same maximum), keep the longest ones.
    subsumed_indices = []
    if len(subsumed_ranges) > 1:
        logger.debug("SUBSUMING")
        for i, r1 in enumerate(subsumed_ranges):
            for j, r2 in enumerate(subsumed_ranges):
                if i != j:
                    if tuple_subsume(r1, r2):
                        subsumed_indices.append(j)
        subsumed_set = set(subsumed_indices)
        kept_indices = [x for x in range(len(subsumed_ranges)) if x not in subsumed_set]
        return kept_indices
    else:
        collapsed_ranges = subsumed_ranges

    if False:  # num_breaks == 0:  # CJJ This isn't used
        kept_indicies = [range_list.index(i) for i in connected_ranges]
        return kept_indicies
    else:
        # List contains other lists, need to flatten this to just tuples.
        flattened_list = []
        for a in range(len(collapsed_ranges)):
            if isinstance(collapsed_ranges[a], list):
                for i in collapsed_ranges[a]:
                    flattened_list.append(i)
            else:
                flattened_list.append(collapsed_ranges[a])
        kept_indicies = [range_list.index(i) for i in flattened_list]
        return kept_indicies


def tuple_overlap(a, b):
    """Given two tuples of length two, determine if the ranges overlap"""
    return a[0] < b[0] < a[1] or b[0] < a[0] < b[1]


def tuple_subsume(a, b):
    """Given two tuples of length two, determine if a has a range that includes b"""
    # if b[0] >= a[0] and b[1] <= a[1]:
    if b[0] > a[0] and b[1] < a[1]:  # CJJ is using >= and <= and there are two good contigs with the same protein
        # hit ranges, they both get removed. We don't want this behaviour.
        return True
    else:
        return False


def reciprocal_best_hit(prot, proteinHits):
    """Given a protein dictionary and the dictionary of all protein dictionaries,
        Return the protein dictionary minus any contigs that have higher percentage hits to other proteins."""
    logger = logging.getLogger("pipeline")

    protname = prot["name"]
    kept_indicies = []

    for contig in prot["assemblyHits"]:
        contigname = contig.split(",")[0]
        contig_idx = prot["assemblyHits"].index(contig)
        maxProt = protname
        for otherProt in proteinHits:
            otherprot_contiglist = [x.split(",")[0] for x in proteinHits[otherProt]["assemblyHits"]]
            if proteinHits[otherProt]["name"] != protname:
                if contigname in otherprot_contiglist:
                    full_contigname = [b for b in proteinHits[otherProt]["assemblyHits"] if contigname in b][0]
                    logger.debug("%s %s" % (contig, full_contigname))
                    otherHit_idx = proteinHits[otherProt]["assemblyHits"].index(full_contigname)

                    target_ranges = [sorted((prot["target_begin"][contig_idx], prot["target_end"][contig_idx])), sorted(
                        (proteinHits[otherProt]["target_begin"][otherHit_idx],
                         proteinHits[otherProt]["target_end"][otherHit_idx]))]
                    logger.debug(repr(target_ranges))
                    # Check that the two contig hits have overlapping ranges.
                    if tuple_overlap(target_ranges[0], target_ranges[1]):
                        logger.debug("%s %s" % (
                        repr(prot["percentid"][contig_idx]), repr(proteinHits[otherProt]["percentid"][otherHit_idx])))
                        if prot["percentid"][contig_idx] < proteinHits[otherProt]["percentid"][otherHit_idx]:
                            logger.debug("contig %s is a better hit to %s" % (contigname, otherProt))
                            maxProt = proteinHits[otherProt]["name"]
                    else:
                        logger.debug("ranges did not overlap")
        if maxProt == protname:
            kept_indicies.append(contig_idx)
    return keep_indicies(kept_indicies, prot)


def paralog_test(exonerate_hits, prot, prefix, paralog_warning_min_cutoff):
    """Gives a warning if there are multiple hits of long length to the same protein"""
    logger = logging.getLogger("pipeline")
    protlength = len(prot)
    hitlengths = [abs(int(x.split(",")[2]) - int(x.split(",")[3])) for x in exonerate_hits["assemblyHits"]]
    logger.debug("protein length: {}".format(protlength))
    logger.debug("Hit lengths:")
    logger.debug(hitlengths)
    # longhits = [x > 0.75 * protlength for x in hitlengths]
    longhits = [x > paralog_warning_min_cutoff * protlength for x in hitlengths]  # CJJ
    if sum(longhits) > 1:
        sys.stderr.write("WARNING: Multiple long-length exonerate hits for {}. Check for paralogs!\n".format(prot.id))
        with open("{}/paralog_warning.txt".format(prefix), 'w') as pw:
            for hit in range(len(exonerate_hits["assemblyHits"])):
                if longhits[hit]:
                    pw.write(prot.id + "\t" + exonerate_hits["assemblyHits"][hit] + "\n")


def myTranslate(nucl):
    """Given a raw sequence of nucleotides, return raw sequence of amino acids."""
    nucseq = Seq(nucl)
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


########################################################################################################################
# Define main(), including argparse options
########################################################################################################################

def main():
    parser = argparse.ArgumentParser(
        description="exonerate_hits.py; Generate gene-by-gene protein and nucleotide files from Bait Capture Assembly")
    parser.add_argument("--debug", help="Print debugging information for development testing.",
                        action="store_true", dest="loglevel", default=False)
    parser.add_argument("proteinfile", help="FASTA file containing one 'bait' sequence per protein.")
    parser.add_argument("assemblyfile", help="FASTA file containing DNA sequence assembly.")
    parser.add_argument("--prefix", help="""Prefix for directory, files, and sequences generated from this assembly. 
            If not specified, will be extracted from assembly file name.""", default=None)
    parser.add_argument("--no_sequences", help="Do not generate protein and nucleotide sequence files.",
                        action="store_true", default=False)
    parser.add_argument("--first_search_filename",
                        help="Location of previously completed Exonerate results. Useful for testing.", default="no")
    parser.add_argument("-t", "--threshold",
                        help="Threshold for Percent Identity between contigs and proteins. default = 55%%", default=55,
                        type=int)
    parser.add_argument("--length_pct",
                        help="Include an exonerate hit if it is at least as long as X percentage of the reference "
                             "protein length. Default = 100%%", default=90, type=int)
    parser.add_argument("--depth_multiplier",
                        help="Accept any full-length hit if it has a coverage depth X times the next best hit. Set to "
                             "zero to not use depth. Default = 10", default=10, type=int)
    parser.add_argument("--nosupercontigs",
                        help="Do not create any supercontigs. The longest single Exonerate hit will be used",
                        action="store_true", dest='nosupercontigs', default=False)  # CJJ
    parser.add_argument("--memory", help="memory (RAM ) to use for bbmap.sh", default=1, type=int)  # CJJ
    parser.add_argument("--threads", help="threads to use for bbmap.sh", default=4, type=int)  # CJJ
    parser.add_argument("--discordant_reads_edit_distance",
                        help="Minimum number of differences between one read of a read pair vs the supercontig "
                             "reference for a read pair to be flagged as discordant", default=7, type=int)  # CJJ
    parser.add_argument("--discordant_reads_cutoff",
                        help="minimum number of discordant reads pairs required to flag a supercontigs as a potential "
                             "hybrid of contigs from multiple paralogs", default=100, type=int)  # CJJ
    parser.add_argument("--paralog_warning_min_cutoff", default=0.75, type=float,
                        help="Minimum length percentage of a contig vs reference protein length for a paralog warning "
                             "to be generated. Default is %(default)s")  # CJJ

    args = parser.parse_args()

    ####################################################################################################################
    # Print run info for script parameters
    ####################################################################################################################
    sys.stderr.write(f'\nCJJ from exonerate_hits.py: {args}\n')
    sys.stderr.write(f'\n memory is : {args.memory} \n')
    sys.stderr.write(f'\n edit distance is : {args.discordant_reads_edit_distance} \n')
    sys.stderr.write(f'\n discordant cutoff is : {args.discordant_reads_cutoff} \n')
    sys.stderr.flush()

    ####################################################################################################################
    # Read in the selected protein reference sequence and SPAdes nucleotide contigs
    ####################################################################################################################
    proteinfilename = args.proteinfile  # CJJ read in protein fasta e.g. At2g47110_baits.fasta
    assemblyfilename = args.assemblyfile  # CJJ read in SPAdes contigs e.g. At2g47110_contigs.fasta

    ####################################################################################################################
    # Create a directory for output files based on the prefix name
    ####################################################################################################################
    if args.prefix:
        prefix = args.prefix
        if os.path.exists(prefix):
            pass
        else:
            os.mkdir(prefix)
    else:
        prefix = os.path.basename(assemblyfilename).split(".")[0]

    ####################################################################################################################
    # If producing supercontigs, identify file of interleaved reads for this sample/gene for mapping
    ####################################################################################################################
    if not args.nosupercontigs:  # CJJ
        # print(f'prefix is: {prefix}')
        gene_folder = os.path.split(prefix)[0]
        # gene_folder = os.path.split(prefix)[1] # FOR TESTING ONLY CHANGE BACK!!!!
        interleaved_reads = f'{gene_folder}/{gene_folder}_interleaved.fasta'
        # interleaved_reads = f'{gene_folder}_interleaved.fasta'  # FOR TESTING ONLY CHANGE BACK!!!!
        # print(interleaved_reads)
    else:
        interleaved_reads = 'None'

    ####################################################################################################################
    # Set up logger
    ####################################################################################################################
    logger = logging.getLogger("pipeline")
    ch = logging.StreamHandler()
    logger.addHandler(ch)
    if args.loglevel:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)

    ####################################################################################################################
    # Read the protein reference and the SPAdes contigs into SeqIO dictionaries
    ####################################################################################################################
    try:
        proteinfile = open(proteinfilename)
    except IOError:
        print("The file %s could not be opened!" % proteinfilename)
        return ()

    try:
        assemblyfile = open(assemblyfilename)
    except IOError:
        print("The file %s could not be opened!" % assemblyfilename)
        return ()
    assembly_dict = SeqIO.to_dict(SeqIO.parse(assemblyfile, 'fasta'))
    protein_dict = SeqIO.to_dict(SeqIO.parse(proteinfile, 'fasta'))

    ####################################################################################################################
    # Run the function initial_exonerate, and sort the SPAdes contig hits
    ####################################################################################################################
    if os.path.exists(args.first_search_filename):  # Shortcut for Testing purposes
        logger.info("Reading initial exonerate results from file {}.".format(first_search_filename))
        sequence_dict = SeqIO.to_dict(SeqIO.parse(first_search_filename, 'fasta'))
    else:
        sequence_dict = initial_exonerate(proteinfilename, assemblyfilename, prefix)
    proteinHits = protein_sort(sequence_dict)
    sys.stderr.write("There were {} exonerate hits for {}.\n".format(len(sequence_dict), proteinfilename))

    ####################################################################################################################
    # Create directories for nucleotide (FNA) and amino acid (FAA) sequences
    ####################################################################################################################
    directory_name = "%s/sequences/FNA" % prefix
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)

    directory_name = "%s/sequences/FAA" % prefix
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)

    ####################################################################################################################
    # Filter the Exonerate SPAdes contig hits
    ####################################################################################################################
    for prot in proteinHits:
        logger.debug(prot)
        logger.debug("Initial hits: %s" % len(proteinHits[prot]["assemblyHits"]))

        ################################################################################################################
        # Perform a paralog test and generate warnings
        ################################################################################################################
        paralog_test(proteinHits[prot], protein_dict[prot], prefix, args.paralog_warning_min_cutoff)

        proteinHits[prot]["reflength"] = len(protein_dict[prot])
        proteinHits[prot] = get_contig_order(proteinHits[prot])
        logger.debug("After get_contig_order: %d" % len(proteinHits[prot]["assemblyHits"]))
        # Remove contigs that are suboptimal hits. Only one protein hit allowed per contig.

        proteinHits[prot] = reciprocal_best_hit(proteinHits[prot], proteinHits)
        #         logger.debug("After RBH: %s" % " ".join(proteinHits[prot]["assemblyHits"]))
        logger.debug("After RBH: %d" % len(proteinHits[prot]["assemblyHits"]))
        if len(proteinHits[prot]["assemblyHits"]) == 0:
            report_no_sequences(proteinHits[prot]["name"])
            continue  # All hits have been filtered out

        # Filter out contigs with a hit below a threshold
        proteinHits[prot] = filter_by_percentid(proteinHits[prot], args.threshold)
        #         logger.debug("After filter_by_percent_id: %s" % " ".join(proteinHits[prot]["assemblyHits"]))
        logger.debug("After filter_by_percent_id: %d" % len(proteinHits[prot]["assemblyHits"]))
        if len(proteinHits[prot]["assemblyHits"]) == 0:
            report_no_sequences(proteinHits[prot]["name"])
            continue  # All hits have been filtered out

        ################################################################################################################
        # Delete contigs if their range is completely subsumed by another hit's range.
        ################################################################################################################
        proteinHits[prot] = overlapping_contigs(proteinHits[prot], args.length_pct * 0.01, args.depth_multiplier)
        logger.debug("After overlapping_contigs: %d" % len(proteinHits[prot]["assemblyHits"]))

        ################################################################################################################
        # Stitch together a "supercontig" containing all the hits and conduct a second exonerate search.
        ################################################################################################################
        if len(proteinHits[prot]["assemblyHits"]) == 0:
            report_no_sequences(proteinHits[prot]["name"])
            continue  # All hits have been filtered out

        nucl_sequence = fullContigs(proteinHits[prot], sequence_dict, assembly_dict, protein_dict, prefix,
                                    args.threshold, args.nosupercontigs, interleaved_reads=interleaved_reads,
                                    memory=args.memory, discordant_cutoff=args.discordant_reads_cutoff,
                                    edit_distance=args.discordant_reads_edit_distance, threads=args.threads)

        ################################################################################################################
        # If a sequence for the locus was returned, translate it, and write nucleotide and protein seqs to file
        ################################################################################################################
        if nucl_sequence:
            if args.no_sequences:
                continue
            else:
                amino_sequence = myTranslate(nucl_sequence)
                seqID = prefix.split("/")[-1].strip("/")
                sys.stderr.write("Writing amino acid sequence, length: {}\n".format(len(amino_sequence)))
                sys.stdout.write("{}\t{}\n".format(prot.split("-")[-1], len(amino_sequence)))
                amino_filename = "%s/sequences/FAA/%s.FAA" % (prefix, prot.split("-")[-1])
                amino_file = open(amino_filename, 'w')
                amino_file.write(">%s\n%s\n" % (seqID, amino_sequence))
                amino_file.close()

                nucleo_filename = "%s/sequences/FNA/%s.FNA" % (prefix, prot.split("-")[-1])
                nucleo_file = open(nucleo_filename, 'w')
                nucleo_file.write(">%s\n%s\n" % (seqID, nucl_sequence))
                nucleo_file.close()
    proteinfile.close()
    assemblyfile.close()

    sys.stderr.write(f"\nCJJ Finishing exonerate_hits.py \n")
    sys.stderr.flush()


########################################################################################################################
# Run the script
########################################################################################################################
if __name__ == "__main__":
    main()

################################################## END OF SCRIPT #######################################################

