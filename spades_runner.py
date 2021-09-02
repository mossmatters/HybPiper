#!/usr/bin/env python

import argparse, os, sys, shutil, subprocess, re

helptext = '''Run the assembler SPAdes with re-dos if any of the k-mers are unsuccessful.
The re-runs are attempted by removing the largest k-mer and re-running spades. If a final
contigs.fasta file is generated, a 'spades.ok' file is saved.'''


def file_exists_and_not_empty(file_name):  # CJJ
    """
    Check if file exists and is not empty by confirming that its size is not 0 bytes
    """
    # Check if file exist and is not empty
    return os.path.isfile(file_name) and not os.path.getsize(file_name) == 0


def make_spades_cmd(genelist, cov_cutoff=8, cpu=None, paired=True, kvals=None, redo=False, timeout=None,
                    unpaired=False, merged=False):
    if kvals:
        kvals = ",".join(kvals)

    parallel_cmd_list = ["time", "parallel", "--eta"]
    if cpu:
        parallel_cmd_list.append("-j {}".format(cpu))
    if timeout:
        parallel_cmd_list.append("--timeout {}%".format(timeout))

    spades_cmd_list = ["spades.py --only-assembler --threads 1 --cov-cutoff", str(cov_cutoff)]
    # spades_cmd_list = ["spades.py --only-assembler --sc --threads 1 --cov-cutoff", str(cov_cutoff)]  # CJJ added --sc
    if kvals:
        spades_cmd_list.append("-k {}".format(kvals))
    if unpaired:
        spades_cmd_list.append("-s {}/{}_unpaired.fasta")
    if paired and not merged:
        spades_cmd_list.append("--12 {}/{}_interleaved.fasta")
    elif not merged:
        spades_cmd_list.append("-s {}/{}_unpaired.fasta")
    if merged:
        spades_cmd_list.append("--merged {}/{}_merged.fastq")
        spades_cmd_list.append("--12 {}/{}_unmerged.fastq")

    # spades_cmd_list.append("-o {{}}/{{}}_spades :::: {} > spades.log".format(genelist))

    # CJJ WRite separate gene name files for those that have merged reads and those that don't. Format the SPAdes
    # command accordingly and run as two separate subprocess commands in spades_initial().
    if merged:
        with open(genelist, 'r') as gene_file:
            contents = gene_file.readlines()
        genelist_list = [gene.strip() for gene in contents]

        genes_with_merged_reads = [gene for gene in genelist_list if file_exists_and_not_empty(
            f'{gene}/{gene}_merged.fastq')]
        with open(f'spades_genelist_with_merged.txt', 'w') as with_merged:
            for gene in genes_with_merged_reads:
                with_merged.write(f'{gene}\n')

        genes_without_merged_reads = [gene for gene in genelist_list if not file_exists_and_not_empty(
            f'{gene}/{gene}_merged.fastq')]
        with open(f'spades_genelist_without_merged.txt', 'w') as without_merged:
            for gene in genes_without_merged_reads:
                without_merged.write(f'{gene}\n')

        sys.stderr.write(f'With merged: {len(genes_with_merged_reads)}\n')
        sys.stderr.flush()

        sys.stderr.write(f'Without merged: {len(genes_without_merged_reads)}\n')
        sys.stderr.flush()

    if merged:
        spades_cmd_list_with_merged = spades_cmd_list.copy()
        spades_cmd_list_with_merged.append("-o {{}}/{{}}_spades :::: {} >> spades.log".format(
            'spades_genelist_with_merged.txt'))
        # print(spades_cmd_list_with_merged)
        spades_cmd_list_without_merged = spades_cmd_list.copy()
        spades_cmd_list_without_merged = [re.sub('--merged {}/{}_merged.fastq', '', item) for item in
                                          spades_cmd_list_without_merged]
        spades_cmd_list_without_merged.append(
            "-o {{}}/{{}}_spades :::: {} >> spades.log".format('spades_genelist_without_merged.txt'))
        # print(spades_cmd_list_without_merged)

        spades_cmd_with_merged = " ".join(parallel_cmd_list) + " " + " ".join(spades_cmd_list_with_merged)
        spades_cmd_without_merged = " ".join(parallel_cmd_list) + " " + " ".join(spades_cmd_list_without_merged)
        # print(spades_cmd_with_merged)
        # print(spades_cmd_without_merged)
        return spades_cmd_with_merged, spades_cmd_without_merged
    else:
        spades_cmd_list.append("-o {{}}/{{}}_spades :::: {} > spades.log".format(genelist))
        spades_cmd = " ".join(parallel_cmd_list) + " " + " ".join(spades_cmd_list)
        return spades_cmd
    # return spades_cmd


def spades_initial(genelist, cov_cutoff=8, cpu=None, paired=True, kvals=None, timeout=None, unpaired=False,
                   merged=False):
    "Run SPAdes on each gene separately using GNU paralell."""
    if os.path.isfile("spades.log"):
        os.remove("spades.log")

    genes = [x.rstrip() for x in open(genelist)]
    # print paired

    if merged:
        spades_cmd_with_merged, spades_cmd_without_merged = make_spades_cmd(
            genelist, cov_cutoff, cpu, paired=paired, kvals=kvals, unpaired=unpaired, merged=merged)
        sys.stderr.write(spades_cmd_with_merged + "\n")
        exitcode = subprocess.call(spades_cmd_with_merged, shell=True)
        if exitcode:
            sys.stderr.write(
                "ERROR: One or more genes had an error with SPAdes assembly. This may be due to low coverage. No "
                "contigs found for the following genes:\n")
        sys.stderr.write(spades_cmd_without_merged + "\n")
        exitcode = subprocess.call(spades_cmd_without_merged, shell=True)
        if exitcode:
            sys.stderr.write(
                "ERROR: One or more genes had an error with SPAdes assembly. This may be due to low coverage. No "
                "contigs found for the following genes:\n")
    else:
        spades_cmd = make_spades_cmd(genelist, cov_cutoff, cpu, paired=paired, kvals=kvals, unpaired=unpaired,
                                     merged=merged)

        sys.stderr.write("Running SPAdes on {} genes\n".format(len(genes)))
        sys.stderr.write(spades_cmd + "\n")
        exitcode = subprocess.call(spades_cmd, shell=True)

        if exitcode:
            sys.stderr.write(
                "ERROR: One or more genes had an error with SPAdes assembly. This may be due to low coverage. No "
                "contigs found for the following genes:\n")

    spades_successful = []
    spades_failed = []

    for gene in genes:
        gene_failed = False
        if os.path.isfile("{}/{}_spades/contigs.fasta".format(gene, gene)):
            contig_file_size = os.stat("{}/{}_spades/contigs.fasta".format(gene, gene)).st_size
            if contig_file_size > 0:
                shutil.copy("{}/{}_spades/contigs.fasta".format(gene, gene), "{}/{}_contigs.fasta".format(gene, gene))
                spades_successful.append(gene)
            else:
                gene_failed = True
        else:
            gene_failed = True

        if gene_failed:
            sys.stderr.write("{}\n".format(gene))
            spades_failed.append(gene)
    return spades_failed


def rerun_spades(genelist, cov_cutoff=8, cpu=None, paired=True):
    genes = [x.rstrip() for x in open(genelist)]

    redo_cmds_file = open("redo_spades_commands.txt", 'w')

    spades_successful = []
    spades_failed = []
    spades_duds = []

    genes_redos = []

    all_redo_kmers = []
    restart_ks = []
    for gene in genes:
        all_kmers = [int(x[1:]) for x in os.listdir(os.path.join(gene, "{}_spades".format(gene))) if x.startswith("K")]
        all_kmers.sort()

        if len(all_kmers) < 2:
            sys.stderr.write("WARNING: All Kmers failed for {}!\n".format(gene))
            spades_duds.append(gene)
            continue
        else:
            genes_redos.append(gene)
        redo_kmers = [str(x) for x in all_kmers[:-1]]
        restart_k = "k{}".format(redo_kmers[-1])
        kvals = ",".join(redo_kmers)
        spades_cmd = "spades.py --restart-from {} -k {} --cov-cutoff {} -o {}/{}_spades".format(restart_k, kvals,
                                                                                                cov_cutoff, gene, gene)
        redo_cmds_file.write(spades_cmd + "\n")

    redo_cmds_file.close()
    if cpu:
        redo_spades_cmd = "parallel -j {} --eta --timeout 400% :::: redo_spades_commands.txt > spades_redo.log".format(
            cpu)
    else:
        redo_spades_cmd = "parallel --eta --timeout 400% :::: redo_spades_commands.txt > spades_redo.log"

    sys.stderr.write("Re-running SPAdes for {} genes\n".format(len(genes_redos)))
    sys.stderr.write(redo_spades_cmd + "\n")
    exitcode = subprocess.call(redo_spades_cmd, shell=True)

    if exitcode:
        sys.stderr.write(
            "ERROR: One or more genes had an error with SPAdes assembly. This may be due to low coverage. No contigs found for the following genes:\n")

    for gene in genes_redos:
        gene_failed = False
        if os.path.isfile("{}/{}_spades/contigs.fasta".format(gene, gene)):
            if os.stat("{}/{}_spades/contigs.fasta".format(gene, gene)).st_size > 0:
                shutil.copy("{}/{}_spades/contigs.fasta".format(gene, gene), "{}/{}_contigs.fasta".format(gene, gene))
                spades_successful.append(gene)
            else:
                gene_failed = True
        else:
            gene_failed = True

        if gene_failed:
            sys.stderr.write("{}\n".format(gene))
            spades_duds.append(gene)
    with open("spades_duds.txt", 'w') as spades_duds_file:
        spades_duds_file.write("\n".join(spades_duds))

    return spades_failed, spades_duds


def main():
    parser = argparse.ArgumentParser(description=helptext, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('genelist',
                        help="Text file containing the name of each gene to conduct SPAdes assembly. One gene per line, should correspond to directories within the current directory.")
    parser.add_argument('--cpu', type=int, default=0,
                        help="Limit the number of CPUs. Default is to use all cores available.")
    parser.add_argument('--cov_cutoff', type=int, default=8, help="Coverage cutoff for SPAdes. default: %(default)s")
    parser.add_argument("--kvals", nargs='+',
                        help="Values of k for SPAdes assemblies. Default is to use SPAdes auto detection based on read lengths (recommended).",
                        default=None)
    parser.add_argument("--redos_only", action="store_true", default=False,
                        help="Continue from previously assembled SPAdes assemblies and only conduct redos from failed_spades.txt")
    parser.add_argument("--single", help="Reads are single end. Default is paired end.", action='store_true',
                        default=False)
    parser.add_argument("--timeout",
                        help="Use GNU Parallel to kill processes that take longer than X times the average.", default=0)
    parser.add_argument("--unpaired", help="For assembly with both paired (interleaved) and unpaired reads",
                        action="store_true", default=False)
    parser.add_argument("--merged", help="For assembly with both merged and unmerged (interleaved) reads",
                        action="store_true", default=False)
    args = parser.parse_args()

    if args.merged:
        sys.stderr.write(f'\n****CJJ MERGE FLAG IS {args.merged} ******************\n')
        sys.stderr.flush()

    if args.unpaired:
        sys.stderr.write(f'\n****CJJ UNPAIRED FLAG IS {args.unpaired} ******************\n')
        sys.stderr.flush()

    if args.single:
        is_paired = False
    else:
        is_paired = True

    if os.path.isfile("failed_spades.txt") and args.redos_only:
        spades_failed = rerun_spades("failed_spades.txt", cpu=args.cpu, paired=is_paired)
    else:
        if args.unpaired:  # Create empty unpaired file if it doesn't exist
            for gene in open(args.genelist):
                gene = gene.rstrip()
                if os.path.isfile("{}/{}_interleaved.fasta".format(gene, gene)):
                    if not os.path.isfile("{}/{}_unpaired.fasta".format(gene, gene)):
                        open("{}/{}_unpaired.fasta".format(gene, gene), 'a').close()

        spades_failed = spades_initial(args.genelist, cov_cutoff=args.cov_cutoff, cpu=args.cpu, kvals=args.kvals,
                                       paired=is_paired, timeout=args.timeout, unpaired=args.unpaired,
                                       merged=args.merged)

        if len(spades_failed) > 0:
            with open("failed_spades.txt", 'w') as failed_spadefile:
                failed_spadefile.write("\n".join(spades_failed))

            spades_failed, spades_duds = rerun_spades("failed_spades.txt", cov_cutoff=args.cov_cutoff, paired=is_paired,
                                                      cpu=args.cpu)
            if len(spades_failed) == 0:
                sys.stderr.write("All redos completed successfully!\n")
            else:
                sys.exit(1)


if __name__ == "__main__": main()
