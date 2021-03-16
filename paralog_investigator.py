#!/usr/bin/env python

helptext = '''This script will extract exonerate results for all competing full-length paralogs
 and deposit them in a new paralog directory within the results directory.
 Paralogs can then be collected using paralog_retriever.py in order to align
 and build gene family trees.'''

import os, sys, argparse
from Bio import SeqIO


def extract_paralogs(gene, prefix):
    putative_paralog_ids = list(
        set([x.split()[1].rstrip() for x in open(os.path.join(gene, prefix, "paralog_warning.txt"))]))
    try:
        chosen_paralog = open(os.path.join(gene, prefix, "exonerate_stats.csv")).readline().rstrip()
    except IOError:
        return 0

    exonerate_dict = SeqIO.to_dict(SeqIO.parse(os.path.join(gene, prefix, "exonerate_results.fasta"), 'fasta'))

    if not os.path.isdir(os.path.join(gene, prefix, 'paralogs')):
        os.mkdir(os.path.join(gene, prefix, "paralogs"))
    seqs_to_write = [exonerate_dict[x] for x in putative_paralog_ids]

    for seq in range(len(seqs_to_write)):
        if seqs_to_write[seq].id == chosen_paralog:
            seqs_to_write[seq].id = "{}.{}".format(prefix, "main")

        else:
            seqs_to_write[seq].id = "{}.{}".format(prefix, seq)

    SeqIO.write(seqs_to_write, os.path.join(gene, prefix, 'paralogs', '{}_paralogs.fasta'.format(gene)), 'fasta')

    return len(seqs_to_write)


def main():
    parser = argparse.ArgumentParser(description=helptext, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('prefix', help='Name of directory containing HybPiper output.')
    parser.add_argument('--genelist',
                        help="Text file containing the name of each gene extract. The default is to use all genes in 'prefix/genes_with_paralog_warnings.txt'.",
                        default=None)

    args = parser.parse_args()

    if args.genelist:
        if os.path.isfile(args.genelist):
            genelist = [x.rstrip() for x in open(args.genelist)]
        else:
            sys.stderr.write("ERROR: cannot find genelist {}\n".format(args.genelist))
            sys.exit(1)
    else:
        if os.path.isfile("{}/genes_with_paralog_warnings.txt".format(args.prefix)):
            genelist = [x.rstrip() for x in open("{}/genes_with_paralog_warnings.txt".format(args.prefix))]
        else:
            sys.stderr.write("ERROR: cannot find genes_with_paralog_warnings.txt in {}!\n".format(args.prefix))
            # sys.exit(1)
            sys.exit(0)

    os.chdir(args.prefix)
    # Now we only need the last component of the prefix path
    prefixParentDir, prefix = os.path.split(args.prefix)
    if not prefix:
        # if prefix has a trailing /, prefixParentDir will have the / stripped and prefix will be empty.
        # so try again
        prefix = os.path.split(prefixParentDir)[1]

    for gene in genelist:
        num_paralogs = extract_paralogs(gene, prefix)
        sys.stderr.write("{} paralogs written for {}\n".format(num_paralogs, gene))


if __name__ == "__main__": main()
