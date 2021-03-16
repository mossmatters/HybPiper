#!/usr/bin/env python

"""
This script will take the output of a run of HybSeqPipeline (exon sequences) and attempt
to extract intron sequences from the SPAdes assemblies. It is important that the
directory structure of HybSeqPipeline was not disturbed, so that it can be used to collect
information and re-execute exonerate.
"""

import sys, os, argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from exonerate_hits import range_connectivity, tuple_overlap


def get_contig_info():
    """
    Given the prefix of a run of HybSeqPipeline, retrieve the stats from Exonerate
    """
    statspath = "exonerate_stats.csv"
    contig_info = [x.rstrip().split(',') for x in open(statspath).readlines()]
    # print contig_info
    sorted_contig_stats = sorted(contig_info, key=lambda x: int(x[2]))
    # print sorted_contig_stats
    return sorted_contig_stats


def make_intron_supercontig(contig_info, gene, prefix, add_N=False):
    """
    CJJ: reads the SPAdes contigs to a dictionary. If a contig had Exonerate hits after filtering (list recovered via
    function get_contig_info() above), concatenate such contigs to create an 'intron supercontig' and write it to a
    fasta file.

    add_N=False by default.
    """
    spades_contigs = SeqIO.to_dict(SeqIO.parse("../{}_contigs.fasta".format(gene), 'fasta'))
    intron_supercontig = SeqRecord(Seq(''))
    for i in contig_info:
        if i[5] == "(+)":
            intron_supercontig += spades_contigs[i[0]]
        elif i[5] == "(-)":
            intron_supercontig += spades_contigs[i[0]].reverse_complement()
        else:
            sys.stderr.write("Strandedness not found!")
            sys.exit(1)
        if add_N and i != contig_info[-1]:
            intron_supercontig += "NNNNNNNNNN"
    intron_supercontig.id = '{}-{}'.format(prefix, gene)
    intron_supercontig.description = ''
    SeqIO.write(intron_supercontig, 'sequences/intron/{}_supercontig.fasta'.format(gene), 'fasta')


def re_run_exonerate(gene, target="new_faa"):
    """
    CJJ: uses the protien gene.FAA sequence created by reads_first.py, and uses it as a query against the
    'intron-supercontig' in a protein2genome Exonerate search. Writes a gff file of results called
    'intronerate_raw.gff'.

    target="new_faa" by default.
    """
    if target == "new_faa":
        exonerate_cmd = "exonerate -m protein2genome -q sequences/FAA/{}.FAA -t sequences/intron/{}_supercontig.fasta " \
                        "--verbose 0 --showalignment no --showvulgar no --showtargetgff yes > " \
                        "intronerate_raw.gff".format(gene, gene)
    else:
        exonerate_cmd = "exonerate -m protein2genome -q ../{}_baits.fasta -t sequences/intron/{}_supercontig.fasta " \
                        "--verbose 0 --showalignment no --showvulgar no --showtargetgff yes > " \
                        "intronerate_raw.gff".format(gene, gene)
    sys.stderr.write("[CMD] {}\n".format(exonerate_cmd))
    os.system(exonerate_cmd)


def parse_gff(filename):
    """
    Parse a GFF file created for a single gene, return a list of lists containing the annotation info.
    """
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


def longest_hit(hits):
    """
    Given a list of hits, return the longest one.

    CJJ: looks like it returns an integer corresponding to an index, to me - is this a bug? Looks like it'll just
    return the last value for 'hrange' - meant to return 'longest_hit'?
    """
    print("Using longest hit for {}\n".format(hits[0][0][0]))
    ranges = [(int(hit[0][3]), int(hit[0][4])) for hit in hits]
    max_length = 0
    for hrange in range(len(ranges)):
        hit_length = ranges[hrange][1] - ranges[hrange][0]
        if hit_length > max_length:
            max_length = hit_length
            longest_hit = hrange
    # return hrange
    return longest_hit  # CJJ


def score_filter(hits, score_multiplier=2):
    """
    Given the GFF hits with overlapping ranges, determine if one has a score far
    exceeding the others and return that one, else return None.
    """
    print("Searching for hit with score {} times better\n".format(score_multiplier))
    scores = [int(x[0][5]) for x in hits]
    # print(scores)
    max_score = max(scores)
    # print(scores.index(max_score))
    for s in scores:
        if s != max_score:
            if s * score_multiplier > max_score:
                print("No top score found")
                return None
    return scores.index(max_score)


def join_zones(hits):
    """
    Join the hits together.
    """
    min_start = 10000000
    max_end = 0
    for h in hits:
        start = int(h[3])
        end = int(h[4])
        if start < min_start:
            min_start = start
        if end > max_end:
            max_end = end
    hits[0][3] = str(min_start)
    hits[0][4] = str(max_end)
    return hits[0]


def filter_gff(hits, merge=True):
    """
    CJJ:
    """
    hits = sorted(hits, key=lambda x: int(x[0][3]))
    # Get only the features annotated as genes
    gene_annotations = [x for y in hits for x in y if x[2] == 'gene']
    # Get the start,end, and score for each gene annotation
    range_list = [(int(x[3]), int(x[4])) for x in gene_annotations]
    # print(range_list)
    kept_indicies = range_connectivity(range_list)
    kept_range_list = [range_list[x] for x in kept_indicies]
    # print(kept_indicies)
    if len(kept_indicies) > 1:
        overlapping_indicies = []
        non_overlapping_indicies = []
        for ix in range(len(kept_indicies) - 1):
            # print range_list[ix],range_list[ix+1],tuple_overlap(range_list[ix],range_list[ix+1])
            if tuple_overlap(kept_range_list[ix], kept_range_list[ix + 1]):
                if kept_indicies[ix] not in overlapping_indicies:
                    #                    if not tuple_overlap(kept_range_list[ix-1],kept_range_list[ix]):
                    overlapping_indicies.append(kept_indicies[ix])
                if kept_indicies[ix + 1] not in overlapping_indicies:
                    overlapping_indicies.append(kept_indicies[ix + 1])
            else:
                non_overlapping_indicies.append(kept_indicies[ix])
        # print overlapping_indicies
        if overlapping_indicies:
            best_score = score_filter([hits[x] for x in overlapping_indicies])
            if best_score:
                non_overlapping_indicies.append(best_score)
            else:
                if merge:
                    # merge the gene first
                    gene_anno = []
                    cds_anno = []
                    exon_anno = []
                    intron_anno = []
                    similarity_anno = []
                    misc_anno = []
                    # print(hits)
                    for x in hits:
                        for l in x:
                            if l[2] == 'gene':
                                gene_anno.append(l)
                            elif l[2] == "cds":
                                cds_anno.append(l)
                            elif l[2] == "exon":
                                exon_anno.append(l)
                            elif l[2] == "intron":
                                intron_anno.append(l)
                            elif l[2] == "similarity":
                                similarity_anno.append(l)
                            else:
                                misc_anno.append(l)
                    print("Merging {} annotations".format(len(gene_anno)))
                    joined_gene = join_zones(gene_anno)
                    joined_cds = join_zones(cds_anno)
                    joined_exon = join_zones(exon_anno)
                    joined_similarity = join_zones(similarity_anno)
                    joined_hit = [joined_gene, joined_cds, joined_exon]
                    if intron_anno:
                        joined_intron = join_zones(intron_anno)
                        joined_hit.append(joined_intron)
                    joined_hit.append(joined_similarity)
                    joined_hit += misc_anno
                    hits.append(joined_hit)
                    kept_indicies.append(len(hits) - 1)
                    non_overlapping_indicies.append(len(kept_indicies) - 1)
                else:
                    longest = longest_hit([hits[x] for x in overlapping_indicies])
                    if longest:
                        non_overlapping_indicies.append(longest)

        #         for pair in overlapping_indicies:
        #             if int(gene_annotations[pair[0]][5]) > int(gene_annotations[pair[1]][5]):
        #                 non_overlapping_indicies.append(pair[0])
        #             else:
        #                 non_overlapping_indicies.append(pair[1])
        #         if not tuple_overlap(range_list[-2],range_list[-1]):
        #                 non_overlapping_indicies.append(kept_indicies[-1])
        # print kept_indicies
        # print non_overlapping_indicies

        return [hits[kept_indicies[x]] for x in sorted(non_overlapping_indicies)]  # .sort()]

    else:
        return [hits[x] for x in kept_indicies]
        # print [gene_annotations[x] for x in kept_indicies]


def get_new_gff(kept_hits):
    """
    CJJ:
    """
    flatter_list = []
    for hit in kept_hits:
        for line in hit:
            flatter_list.append("\t".join(line))
    return "\n".join(flatter_list) + '\n'


def remove_exons(gff_filename, supercontig_filename, mode="all"):
    """
    Given a supercontig and corresponding annotation, remove the exon sequences. In "intron" mode, only return
    sequences specifically annotated as introns.
    """
    exon_starts = []
    exon_ends = []
    gff = open(gff_filename).readlines()
    for line in gff:
        line = line.rstrip().split("\t")
        if len(line) > 2:
            if line[2] == "exon":
                exon_starts.append(int(line[3]))
                exon_ends.append(int(line[4]))
    supercontig = SeqIO.read(supercontig_filename, 'fasta')
    exonless_contig = SeqRecord(Seq(''), id=supercontig.id)
    start = 0
    for exon in range(len(exon_starts)):
        exonless_contig += supercontig[start:exon_starts[exon] - 1]
        start = exon_ends[exon]
    exonless_contig += supercontig[start:]
    exonless_contig.description = ''
    return exonless_contig


def check_for_files(gene, prefix):
    """
    Check to see if the files needed for intronerate are really present.
    """
    if os.path.isfile("{}/{}/exonerate_stats.csv".format(gene, prefix)):
        if os.path.isfile("{}/{}/sequences/FAA/{}.FAA".format(gene, prefix, gene)):
            if os.path.isfile("{}/{}_contigs.fasta".format(gene, gene)):
                return True
            else:
                return False
        else:
            return False
    else:
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--genelist",
                        help="Optional list of genes to retrieve coverage. Default is to use genes_with_seqs.txt")
    parser.add_argument("--prefix", help="Prefix of sample directory generated by HybSeqPipeline", required=True)
    parser.add_argument("--no-exonerate", help="Don't re-run exonerate, use existing intronerate gff files.",
                        action='store_true', default=False)
    parser.add_argument("--use_target",
                        help="Align the supercontig to the original target sequences, rather than the newly generated "
                             "FAA", default=False, action="store_true")
    parser.add_argument("--merge",
                        help="Merge overlapping annotations for genes, exons, and introns. Default is to pick the "
                             "longest annotation", action="store_true", default=False)
    parser.add_argument("--addN",
                        help="Insert 10 Ns between the contigs when constructing the supercontig. Useful to identify "
                             "where intron recovery fails.", default=False, action="store_true")
    # parser.add_argument("--introns-only",help = "In the intron.fasta file for each gene, only write regions annotated
    # as introns by exonerate. Default: all non-exon regions are written to introns.fasta.",action="store_true",
    # default=False)

    args = parser.parse_args()

    if args.genelist:  # CJJ: defaults to file 'genes_with_seqs.txt' created by reads_first.py otherwise.
        genelist_fn = os.path.abspath(args.genelist)

    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    if os.path.isdir(args.prefix):
        os.chdir(args.prefix)
        basedir = os.getcwd()
        prefix = os.path.split(basedir)[1]
    else:
        sys.stderr.write("Directory {} not found!\n".format(args.prefix))

    ####################################################################################################################
    # CJJ: recover a list of gene names that have seqs recovered via reads_first.py (listed in file
    # 'genes_with_seqs.txt')
    ####################################################################################################################
    if args.genelist:
        genelist = [x.split()[0] for x in open(genelist_fn).readlines()]
    else:
        try:
            genelist = [x.split()[0] for x in open('genes_with_seqs.txt').readlines()]
        except FileNotFoundError:
            sys.stderr.write(f'No "genes_with_seqs.txt" file found, exiting!')
            sys.exit(0)

    ####################################################################################################################
    # CJJ: placeholder
    ####################################################################################################################
    with open("intron_stats.txt", 'w') as intron_stats_file:
        full_gff = ''
        for gene in genelist:
            if check_for_files(gene, prefix):
                os.chdir("{}/{}".format(gene, prefix))
                contig_info = get_contig_info()  # CJJ: get a list of SPAdes contigs that had Exonerate hits after
                # filtering.
                if not os.path.exists("sequences/intron"):
                    os.makedirs("sequences/intron")
                make_intron_supercontig(contig_info, gene, prefix, add_N=args.addN)  # CJJ: write a fasta file of
                # concatenated SPADes contigs (i.e. those that had Exonerate hits after filtering).
                if not args.no_exonerate:
                    if args.use_target:
                        re_run_exonerate(gene, target="bait")
                    else:
                        re_run_exonerate(gene, target="new_faa")  # CJJ: default via argparse
                hits = parse_gff("intronerate_raw.gff")  # CJJ: from a gff Exonerate result, return a list of lists
                # containing the annotation info (split on tabs?).

                # print([h[0] for h in hits])
                kept_hits = filter_gff(hits, merge=args.merge)  # CJJ: args.merge is False by default. Default is to
                # pick the longest annotation.

                # print([h[0] for h in kept_hits])
                with open("intronerate.gff", 'w') as new_gff:
                    new_gff_string = get_new_gff(kept_hits)
                    num_introns = new_gff_string.count("intron\t")
                    intron_stats_file.write("{}\t{}\t{}\n".format(prefix, gene, num_introns))
                    sys.stderr.write("{} introns found for {}.\n".format(num_introns, gene))

                    new_gff.write(new_gff_string)
                    full_gff += new_gff_string
                exonless_contig = remove_exons("intronerate.gff", "sequences/intron/{}_supercontig.fasta".format(gene))
                SeqIO.write(exonless_contig, "sequences/intron/{}_introns.fasta".format(gene), 'fasta')

                os.chdir(basedir)
            else:
                sys.stderr.write("ERROR: Files not found for {} gene {}!\n".format(prefix, gene))
    with open("{}_genes.gff".format(prefix), 'w') as all_gff:
        all_gff.write(full_gff)


if __name__ == "__main__":
    main()
