#!/usr/bin/env python

import argparse,os,sys,importlib,shutil,subprocess,glob

helptext="""

HybPiper Version 1.2 (March 2017)

This script is a wrapper around several scripts in the HybSeqPipeline.
It can check whether you have the appropriate dependencies available (see --check-depend).
It makes sure that the other scripts needed are in the same directory as this one.
Command line options are passed to the other executables.
Unless --prefix is set, output will be put within a directory named after your read files."""

velvet_genefilename = "velvet_genelist.txt"
cap3_genefilename = "cap3_genelist.txt"
exonerate_genefilename = "exonerate_genelist.txt"
spades_genefilename = "spades_genelist.txt"


def py_which(cmd, mode=os.F_OK | os.X_OK, path=None):
    """Given a command, mode, and a PATH string, return the path which
    conforms to the given mode on the PATH, or None if there is no such
    file.

    `mode` defaults to os.F_OK | os.X_OK. `path` defaults to the result
    of os.environ.get("PATH"), or can be overridden with a custom search
    path.

    """
    # Check that a given file can be accessed with the correct mode.
    # Additionally check that `file` is not a directory, as on Windows
    # directories pass the os.access check.
    def _access_check(fn, mode):
        return (os.path.exists(fn) and os.access(fn, mode)
                and not os.path.isdir(fn))

    # If we're given a path with a directory part, look it up directly rather
    # than referring to PATH directories. This includes checking relative to the
    # current directory, e.g. ./script
    if os.path.dirname(cmd):
        if _access_check(cmd, mode):
            return cmd
        return None

    if path is None:
        path = os.environ.get("PATH", os.defpath)
    if not path:
        return None
    path = path.split(os.pathsep)

    if sys.platform == "win32":
        # The current directory takes precedence on Windows.
        if not os.curdir in path:
            path.insert(0, os.curdir)

        # PATHEXT is necessary to check on Windows.
        pathext = os.environ.get("PATHEXT", "").split(os.pathsep)
        # See if the given file matches any of the expected path extensions.
        # This will allow us to short circuit when given "python.exe".
        # If it does match, only test that one, otherwise we have to try
        # others.
        if any([cmd.lower().endswith(ext.lower()) for ext in pathext]):
            files = [cmd]
        else:
            files = [cmd + ext for ext in pathext]
    else:
        # On other platforms you don't have things like PATHEXT to tell you
        # what file suffixes are executable, so just pass on cmd as-is.
        files = [cmd]

    seen = set()
    for dir in path:
        normdir = os.path.normcase(dir)
        if not normdir in seen:
            seen.add(normdir)
            for thefile in files:
                name = os.path.join(dir, thefile)
                if _access_check(name, mode):
                    return name
    return None



def check_dependencies():
    """Checks for the presence of executables and Python packages"""
    executables =         ["blastx",
                        "exonerate",
                        #"velvetg",
                        #"velveth",
                        #"cap3",
                        "parallel",
                        "makeblastdb",
                        "spades.py",
                        "bwa",
                        "samtools"]

    python_packages =     ["Bio"]
    
    everything_is_awesome = True
    for e in executables:
        e_loc = py_which(e)
        if e_loc:
            print(("{} found at {}".format(e,e_loc)))
        else:
            print(("{} not found in your $PATH!".format(e)))
            everything_is_awesome = False
            
    for p in python_packages:
        try:
            i = importlib.import_module(p)
            print(("Package {} successfully loaded!".format(p)))
        except ImportError:
            print(("Package {} not found!".format(p)))
            everything_is_awesome=False
    return everything_is_awesome


def blastx(readfiles,baitfile,evalue,basename,cpu=None,max_target_seqs=10,unpaired=False):
    dna = set("ATCGN")
    if os.path.isfile(baitfile):
        #Quick detection of whether baitfile is DNA. 
        with open(baitfile) as bf:
            header = bf.readline()
            seqline = bf.readline().rstrip().upper()
            if not set(seqline) - dna:
                print("ERROR: only ATCGN characters found in first line. You need a protein bait file for BLASTx!")
                return None

        if os.path.isfile(os.path.split(baitfile)[0]+'.psq'):
            db_file = baitfile
        else:
            print("Making protein blastdb in current directory.")
            if os.path.split(baitfile)[0]:
                shutil.copy(baitfile,'.')
            db_file = os.path.split(baitfile)[1]
            makeblastdb_cmd = "makeblastdb -dbtype prot -in {}".format(db_file)
            print(makeblastdb_cmd)
            exitcode = subprocess.call(makeblastdb_cmd,shell=True)
            if exitcode:
                return None
    else:
        print(("Cannot find baitfile at: {}".format(baitfile)))
        return None

    #Remove previous blast results if they exist (because we will be appending)
    if os.path.isfile(basename+".blastx"):
        os.remove(basename+".blastx")
    
    if unpaired:
        read_file = readfiles
        pipe_cmd = "cat {} |  awk '{{if(NR % 4 == 1 || NR % 4 == 2) {{sub(/@/, \">\"); print; }} }}'".format(read_file)
        blastx_command = "blastx -db {} -query - -evalue {} -outfmt 6 -max_target_seqs {}".format(db_file,evalue,max_target_seqs)
        if cpu:
            full_command = "time {} | parallel -j {} -k --block 200K --recstart '>' --pipe '{}' >> {}_unpaired.blastx ".format(pipe_cmd,cpu,blastx_command,basename)
        else:
            full_command = "time {} | parallel -k --block 200K --recstart '>' --pipe '{}' >> {}_unpaired.blastx ".format(pipe_cmd,blastx_command,basename)
        print(full_command)
        exitcode = subprocess.call(full_command,shell=True)
        if exitcode:
            #Concatenate the two blastfiles.
        
            return None
        return basename + "_unpaired.blastx"
        
    else:
        for read_file in readfiles:

        #Piping commands for Fastq -> FASTA    
        # Curly braces must be doubled within a formatted string.
            pipe_cmd = "cat {} |  awk '{{if(NR % 4 == 1 || NR % 4 == 2) {{sub(/@/, \">\"); print; }} }}'".format(read_file)
    
            blastx_command = "blastx -db {} -query - -evalue {} -outfmt 6 -max_target_seqs {}".format(db_file,evalue,max_target_seqs)
            if cpu:
                full_command = "time {} | parallel -j {} -k --block 200K --recstart '>' --pipe '{}' >> {}.blastx ".format(pipe_cmd,cpu,blastx_command,basename)
            else:
                full_command = "time {} | parallel -k --block 200K --recstart '>' --pipe '{}' >> {}.blastx ".format(pipe_cmd,blastx_command,basename)
            print(full_command)
            exitcode = subprocess.call(full_command,shell=True)
            if exitcode:
                #Concatenate the two blastfiles.
        
                return None
            
    return basename + '.blastx'

def distribute(blastx_outputfile,readfiles,baitfile,run_dir,target=None,unpaired_readfile=None,exclude=None):
    #NEED TO ADD SOMETHING ABOUT DIRECTORIES HERE.
    #print run_dir
    read_cmd = "time python {} {} {}".format(os.path.join(run_dir,"distribute_reads_to_targets.py"),blastx_outputfile," ".join(readfiles))
    exitcode = subprocess.call(read_cmd,shell=True)
    if exitcode:
        print("ERROR: Something went wrong with distributing reads to gene directories.")
        return exitcode
    target_cmds = ["time python", os.path.join(run_dir,"distribute_targets.py"),baitfile,"--blastx",blastx_outputfile]
    if target:
        target_cmds.append("--target {}".format(target))
    if exclude:
        target_cmds.append("-- exclude {}".format(exclude))
    if unpaired_readfile:
        blastx_outputfile = blastx_outputfile.replace(".blastx","_unpaired.blastx")
        unpaired_cmd = "time python {} {} {}".format(os.path.join(run_dir,"distribute_reads_to_targets.py"),blastx_outputfile,unpaired_readfile)
        print(("[CMD] {}\n".format(unpaired_cmd)))
        exitcode = subprocess.call(unpaired_cmd,shell=True)
    target_cmd = " ".join(target_cmds)
    exitcode = subprocess.call(target_cmd,shell=True)
    if exitcode:
        print("ERROR: Something went wrong distributing targets to gene directories.")
        return exitcode
    return None

def distribute_bwa(bamfile,readfiles,baitfile,run_dir,target=None,unpaired=None,exclude=None):
    #NEED TO ADD SOMETHING ABOUT DIRECTORIES HERE.
    #print run_dir
    read_cmd = "time python {} {} {}".format(os.path.join(run_dir,"distribute_reads_to_targets_bwa.py"),bamfile," ".join(readfiles))
    print(("[CMD] {}\n".format(read_cmd)))
    exitcode = subprocess.call(read_cmd,shell=True)
    
    if unpaired:
        up_bamfile = bamfile.replace(".bam","_unpaired.bam")
        unpaired_cmd = "time python {} {} {}".format(os.path.join(run_dir,"distribute_reads_to_targets_bwa.py"),up_bamfile,unpaired)
        print(("[CMD] {}\n".format(unpaired_cmd)))
        exitcode = subprocess.call(unpaired_cmd,shell=True)
    
    if exitcode:
        print("ERROR: Something went wrong with distributing reads to gene directories.")
        return exitcode
    target_cmds = ["time python", os.path.join(run_dir,"distribute_targets.py"),baitfile,"--bam",bamfile]
    if target:
        target_cmds.append("--target {}".format(target))
    if unpaired:
        target_cmds.append("--unpaired")
    if exclude:
        target_cmds.append("--exclude {}".format(exclude))
    target_cmd = " ".join(target_cmds)
    print("[DISTRIBUTE]: {}".format(target_cmd))
    exitcode = subprocess.call(target_cmd,shell=True)
    if exitcode:
        print("ERROR: Something went wrong distributing targets to gene directories.")
        return exitcode
    return None



def make_basename(readfiles,prefix=None):
    """Unless prefix is set, generate a directory based off the readfiles, using everything up to the first underscore.
        If prefix is set, generate the directory "prefix" and set basename to be the last component of the path.
        
        """
    if prefix:
            if not os.path.exists(prefix):
                os.makedirs(prefix)
#             else:
#                 if os.listdir(prefix):
#                 ## When running under HTCondor, the execution might be restarted
#                 ## However, this code in not tolerant of restarts
#                     print "ERROR: Directory {} not empty, exiting!".format(prefix)
#                     return

            prefixParentDir, prefix = os.path.split(prefix)
            if not prefix:
                # if prefix has a trailing /, prefixParentDir will have the / stripped and prefix will be empty.
                # so try again
                prefix = os.path.split(prefixParentDir)[1]


            return prefixParentDir,prefix

        ## --prefix is not set on cmd line;  Write output to subdir in . 
    basename = os.path.split(readfiles[0])[1].split('_')[0]

    if not os.path.exists(basename):
        os.makedirs(basename)
    return '.', basename

def spades(genes,run_dir,cov_cutoff=8,cpu=None,paired=True,kvals=None,timeout=None,unpaired=False):
    "Run SPAdes on each gene separately using GNU paralell."""
    
#    import spades_runner

    with open(spades_genefilename,'w') as spadesfile:
        spadesfile.write("\n".join(genes)+"\n")
    
    if os.path.isfile("spades.log"):
        os.remove("spades.log")
    if os.path.isfile("spades_redo.log"):
        os.remove("spades_redo.log")

    spades_runner_list = ["python","{}/spades_runner.py".format(run_dir),spades_genefilename,"--cov_cutoff",str(cov_cutoff)]
    if cpu:
        spades_runner_list.append("--cpu")
        spades_runner_list.append(str(cpu))
    if not paired:
        spades_runner_list.append("--single")
    if unpaired:
        spades_runner_list.append("--unpaired")
    if timeout:
        spades_runner_list.append("--timeout")
        spades_runner_list.append("{}%".format(timeout))
    if kvals:
        spades_runner_list.append("--kvals")
        spades_runner_list.append("{}".format(",".join(kvals)))
    
    
    spades_runner_cmd = " ".join(spades_runner_list)    

#     if cpu:
#         if paired:
#             spades_runner_cmd = "python {} {} --cpu {} --cov_cutoff {}".format(os.path.join(run_dir,"spades_runner.py"),spades_genefilename, cpu, cov_cutoff)
#         else:
#             spades_runner_cmd = "python {} {} --cpu {} --cov_cutoff {} --single".format(os.path.join(run_dir,"spades_runner.py"),spades_genefilename, cpu, cov_cutoff)
#     else:
#         if paired:
#             spades_runner_cmd = "python {} {} --cov_cutoff {}".format(os.path.join(run_dir,"spades_runner.py"),spades_genefilename, cov_cutoff)
#         else:
#             spades_runner_cmd = "python {} {} --cov_cutoff {} --single".format(os.path.join(run_dir,"spades_runner.py"),spades_genefilename, cov_cutoff)
    exitcode = subprocess.call(spades_runner_cmd,shell=True)
    if exitcode:
        sys.stderr.write("WARNING: Something went wrong with the assemblies! Check for failed assemblies and re-run! \n")
        return None
    else:
        if os.path.isfile("spades_duds.txt"):
            spades_duds = [x.rstrip() for x in open("spades_duds.txt")]
        else:
            spades_duds = []



# 
#     spades_failed = spades_runner.spades_initial(spades_genefilename,cov_cutoff,cpu,kvals)
#     if len(spades_failed) > 0:
#         with open("failed_spades.txt",'w') as failed_spadefile:
#             failed_spadefile.write("\n".join(spades_failed))
#     
#         spades_failed_redos,spades_duds = spades_runner.rerun_spades("failed_spades.txt")
#         if len(spades_failed) == 0:
#             sys.stderr.write("All redos completed successfully!\n")
#     
    spades_genelist = []
    for gene in genes:
#        if gene not in set(spades_failed):
        if gene not in set(spades_duds):
#                if gene not in set(spades_failed_redos):
            spades_genelist.append(gene)
    
    with open(exonerate_genefilename,'w') as genefile:
        genefile.write("\n".join(spades_genelist)+"\n")

    return spades_genelist    

# def spades(genes,cov_cutoff=8,cpu=None,paired=True,kvals=None):
#     "Run SPAdes on each gene separately using GNU paralell."""
#     if os.path.isfile("spades.log"):
#         os.remove("spades.log")
#     
#     with open(spades_genefilename,'w') as spadesfile:
#         spadesfile.write("\n".join(genes)+"\n")
#     
#     if paired:
#         fileflag = "--12"
#     else:
#         fileflag = "-s"
#     
#     if kvals:
#         kvals = ",".join(kvals)
#     
#     if cpu:
#         if kvals:
#             spades_cmd = "time parallel -j {} --eta spades.py --only-assembler -k {} --threads 1 --cov-cutoff {} {} {{}}/{{}}_interleaved.fasta -o {{}}/{{}}_spades :::: {} > spades.log".format(cpu,kvals,cov_cutoff,fileflag,spades_genefilename)
#         else:
#             spades_cmd = "time parallel -j {} --eta spades.py --only-assembler --threads 1 --cov-cutoff {} {} {{}}/{{}}_interleaved.fasta -o {{}}/{{}}_spades :::: {} > spades.log".format(cpu,cov_cutoff,fileflag,spades_genefilename)
#     else:
#         if kvals:
#             spades_cmd = "time parallel --eta spades.py --only-assembler -k {} --threads 1 --cov-cutoff {} {} {{}}/{{}}_interleaved.fasta -o {{}}/{{}}_spades :::: {} > spades.log".format(kvals,cov_cutoff,fileflag,spades_genefilename)
#         else:
#             spades_cmd = "time parallel --eta spades.py --only-assembler --threads 1 --cov-cutoff {} {} {{}}/{{}}_interleaved.fasta -o {{}}/{{}}_spades :::: {} > spades.log".format(cov_cutoff,fileflag,spades_genefilename)
#         
#     sys.stderr.write("Running SPAdes on {} genes\n".format(len(genes)))
#     sys.stderr.write(spades_cmd + "\n")
#     exitcode = subprocess.call(spades_cmd,shell=True)
#     
#     #Need to handle an error differently with SPAdes, which can fail if there are simply not enough reads. 
#     
#     if exitcode:
#         sys.stderr.write("ERROR: One or more genes had an error with SPAdes assembly. This may be due to low coverage. No contigs found for the following genes:\n")
#     
#     spades_genelist = []
#     
#     for gene in genes:
#         if os.path.isfile("{}/{}_spades/contigs.fasta".format(gene,gene)):
#             shutil.copy("{}/{}_spades/contigs.fasta".format(gene,gene),"{}/{}_contigs.fasta".format(gene,gene))
#             spades_genelist.append(gene)
#         else:
#             sys.stderr.write("{}\n".format(gene))
#     
#     with open(exonerate_genefilename,'w') as genefile:
#         genefile.write("\n".join(spades_genelist)+"\n")
# 
#     return spades_genelist    
        

def velvet(genes,cov_cutoff=5,ins_length=200,kvals = ["21","31","41","51","61"],cpu=None,paired=True):
    """Use parallel to run velveth and velvetg on a set of k values on every gene with blastx hits from the previous steps."""
    if os.path.isfile('velveth.log'):
        os.remove('velveth.log')
    if os.path.isfile('velvetg.log'):
        os.remove('velvetg.log')
    
    #Write a file with the list of genes for velvet, to avoid screen spam and hitting the argument cap.
    with open(velvet_genefilename,'w') as genefile:
        genefile.write("\n".join(genes)+"\n")
    
    if paired:
        if cpu:
            velveth_cmd = "time parallel  -j {} --eta velveth {{1}}/velvet{{2}} {{2}} -shortPaired {{1}}/{{1}}_interleaved.fasta '>>' velveth.log :::: {} ::: {}".format(cpu,velvet_genefilename," ".join(kvals))
        else:
            velveth_cmd = "time parallel  --eta velveth {{1}}/velvet{{2}} {{2}} -shortPaired {{1}}/{{1}}_interleaved.fasta '>>' velveth.log :::: {} ::: {}".format(velvet_genefilename," ".join(kvals))
    else:
        if cpu:
            velveth_cmd = "time parallel  -j {} --eta velveth {{1}}/velvet{{2}} {{2}} -short {{1}}/{{1}}_interleaved.fasta '>>' velveth.log :::: {} ::: {}".format(cpu,velvet_genefilename," ".join(kvals))
        else:
            velveth_cmd = "time parallel  --eta velveth {{1}}/velvet{{2}} {{2}} -short {{1}}/{{1}}_interleaved.fasta '>>' velveth.log :::: {} ::: {}".format(velvet_genefilename," ".join(kvals))
    

    #print os.getcwd()
    #print os.listdir(".")
    
    
    print(("Running velveth on {} genes".format(len(genes))))
    print(velveth_cmd)
    exitcode = subprocess.call(velveth_cmd,shell=True)
    if exitcode:
        print("ERROR: Something went wrong with velveth!")
        return exitcode
    if cpu:
        velvetg_cmd = "time parallel -j {} --eta velvetg {{1}}/velvet{{2}} -ins_length {} -cov_cutoff {} '>>' velvetg.log :::: {} ::: {}".format(cpu,ins_length,cov_cutoff,velvet_genefilename," ".join(kvals))
    else:
        velvetg_cmd = "time parallel --eta velvetg {{1}}/velvet{{2}} -ins_length {} -cov_cutoff {} '>>' velvetg.log :::: {} ::: {}".format(ins_length,cov_cutoff,velvet_genefilename," ".join(kvals))
    print(("Running velvetg on {} genes".format(len(genes))))
    print(velvetg_cmd)
    exitcode = subprocess.call(velvetg_cmd,shell=True)
    if exitcode:
        print("ERROR: Something went wrong with velvetg!")
        return exitcode
    with open(exonerate_genefilename,'w') as genefile:
        genefile.write("\n".join([x for x in genes if os.path.getsize(os.path.join(x,'velvet_contigs.fa')) > 0]))    
    return None
    
def cap3(genes,cpu=None):
    print("Concatenating velvet output")
    
    if cpu:
        cat_cmd = "time parallel -j {} cat {{1}}/*/contigs.fa '>' {{1}}/velvet_contigs.fa :::: {}".format(cpu,velvet_genefilename)
    else:
        cat_cmd = "time parallel cat {{1}}/*/contigs.fa '>' {{1}}/velvet_contigs.fa :::: {}".format(velvet_genefilename)
    print(cat_cmd)
    exitcode = subprocess.call(cat_cmd,shell=True)
    if exitcode:
        print("ERROR: Something went wrong while concatenating the velvet output files.")
        return exitcode
    #Check that the velvet output actually has data in it.
    genes = [x for x in genes if os.path.getsize(os.path.join(x,'velvet_contigs.fa')) > 0]
    
    with open(cap3_genefilename,'w') as genefile:
        genefile.write("\n".join(genes)+"\n")
    
    if len(genes) == 0:
        print("No genes left! Exiting!")
        return 1
        
    print(("Running Cap3 on {} genes".format(len(genes))))
    if cpu:
        cap3_cmd = "time parallel -j {} --eta cap3 {{1}}/velvet_contigs.fa -o 20 -p 99 '>' {{1}}/{{1}}_cap3.log :::: {}".format(cpu,cap3_genefilename)
    else:
        cap3_cmd = "time parallel --eta cap3 {{1}}/velvet_contigs.fa -o 20 -p 99 '>' {{1}}/{{1}}_cap3.log :::: {}".format(cap3_genefilename)
    print(cap3_cmd)
    exitcode = subprocess.call(cap3_cmd,shell=True)
    if exitcode:
        print("ERROR: Something went wrong with CAP3!")
        return exitcode

    print("Joining CAP3 contigs and singletons")
    join_cmd = "time parallel cat {{1}}/velvet_contigs.fa.cap.contigs {{1}}/velvet_contigs.fa.cap.singlets '>' {{1}}/{{1}}_cap3ed.fa :::: {}".format(cap3_genefilename)
    print(join_cmd)
    exitcode = subprocess.call(join_cmd,shell=True)
    if exitcode:
        print("ERROR: Something went wrong joining the CAP3 output files!")
        return exitcode
        
    return None    

def exonerate(genes,basename,run_dir,replace=True,cpu=None,thresh=55,use_velvet=False,depth_multiplier=0,length_pct=100,timeout=None):
    #Check that each gene in genes actually has CAP3 output
    #cap3_sizes = [os.stat(os.path.join(x,x+"_cap3ed.fa")).st_size for x in genes]
    #print cap3_sizes
    
    if replace:
        for g in genes:
            if os.path.isdir(os.path.join(g,basename)):
                shutil.rmtree(os.path.join(g,basename))
    #genes = [x for x in genes if os.stat(os.path.join(x,x+"_cap3ed.fa")).st_size > 0]
    if len(genes) == 0:
        print(("ERROR: No genes recovered for {}!".format(basename)))
        return 1
    
    if os.path.isfile("genes_with_seqs.txt"):
        os.remove("genes_with_seqs.txt")
    
    if use_velvet:
        file_stem = "cap3ed.fa"
    else:
        file_stem = "contigs.fasta"
    
    print(("Running Exonerate to generate sequences for {} genes".format(len(genes))))
    
    parallel_cmd_list = ["time parallel","--eta"]
    if cpu:
        parallel_cmd_list.append("-j {}".format(cpu))
    if timeout:
        parallel_cmd_list.append("--timeout {}%".format(timeout))    
    
    exonerate_cmd_list = ["python","{}/exonerate_hits.py".format(run_dir),
                "{}/{}_baits.fasta","{{}}/{{}}_{}".format(file_stem),
                "--prefix {{}}/{}".format(basename), 
                "-t {}".format(thresh),
                "--depth_multiplier {}".format(depth_multiplier),
                "--length_pct {}".format(length_pct),
                "::::",
                exonerate_genefilename,
                "> genes_with_seqs.txt"]
                
    exonerate_cmd = " ".join(parallel_cmd_list) + " " + " ".join(exonerate_cmd_list)            
    
#     if cpu:
#         exonerate_cmd = "time parallel -j {} python {} {{}}/{{}}_baits.fasta {{}}/{{}}_{} --prefix {{}}/{} -t {} --depth_multiplier {} --length_pct {} :::: {} > genes_with_seqs.txt".format(cpu,os.path.join(run_dir,"exonerate_hits.py"),file_stem,basename,thresh,depth_multiplier,length_pct,exonerate_genefilename)
#     else:
#         exonerate_cmd = "time parallel python {} {{}}/{{}}_baits.fasta {{}}/{{}}_{} --prefix {{}}/{} -t {} --depth_multiplier {} --length_pct {} :::: {} > genes_with_seqs.txt".format(os.path.join(run_dir,"exonerate_hits.py"),file_stem,basename,thresh,depth_multiplier,length_pct,exonerate_genefilename)
    print(exonerate_cmd)
    exitcode = subprocess.call(exonerate_cmd,shell=True)
    if exitcode:
        print("ERROR: Something went wrong with Exonerate!")
        return exitcode
    return

def bwa(readfiles,baitfile,basename,cpu,unpaired=None):
    """Conduct BWA search of reads against the baitfile.
    Returns an error if the second line of the baitfile contains characters other than ACTGN"""
    dna = set("ATCGN")
    if os.path.isfile(baitfile):
        #Quick detection of whether baitfile is DNA. 
        with open(baitfile) as bf:
            header = bf.readline()
            seqline = bf.readline().rstrip().upper()
            if set(seqline) - dna:
                print("ERROR: characters other than ACTGN found in first line. You need a nucleotide bait file for BWA!")
                return None
    
        if os.path.isfile(os.path.split(baitfile)[0]+'.amb'):
            db_file = baitfile
        else:
            print("Making nucleotide bwa index in current directory.")
            baitfileDir = os.path.split(baitfile)[0]
            if baitfileDir:
                if os.path.realpath(baitfileDir) != os.path.realpath('.'):
                    shutil.copy(baitfile,'.')
            db_file = os.path.split(baitfile)[1]
            make_bwa_index_cmd = "bwa index {}".format(db_file)
            print(("[CMD]: {}".format(make_bwa_index_cmd)))
            exitcode = subprocess.call(make_bwa_index_cmd,shell=True)
            if exitcode:
                return None
    else:
        print(("ERROR: Cannot find baitfile at: {}".format(baitfile)))
        return None
    
    if not cpu:
        import multiprocessing
        cpu = multiprocessing.cpu_count()
    
    if len(readfiles) < 3:
        bwa_fastq = " ".join(readfiles)
    else:
        bwa_fastq = readfiles
            
    bwa_commands = ["time bwa mem","-t",str(cpu),db_file,bwa_fastq," | samtools view -h -b -S - > "]
    if unpaired:
        bwa_commands.append(basename+"_unpaired.bam")
    else:
        bwa_commands.append(basename+".bam")
    full_command = " ".join(bwa_commands)
    print(("[CMD]: {}".format(full_command)))
    exitcode = subprocess.call(full_command,shell=True)
    if exitcode:
        return None
            
    return basename + '.bam'


            
def main():
    parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--check-depend",dest='check_depend',help="Check for dependencies (executables and Python packages) and exit. May not work at all on Windows.",action='store_true')
    parser.add_argument("--bwa",dest="bwa",action='store_true',help="Use BWA to search reads for hits to target. Requires BWA and a bait file that is nucleotides!",default=False)
    parser.add_argument("--no-blast",dest="blast",action="store_false",help="Do not run the blast step. Downstream steps will still depend on the *_all.blastx file. \nUseful for re-runnning assembly/exonerate steps with different options.")
    parser.add_argument("--no-distribute",dest="distribute",action="store_false",help="Do not distribute the reads and bait sequences to sub-directories.")
    parser.add_argument("--no-velvet",dest="velvet",action="store_false",help="Do not run the velvet stages (velveth and velvetg)")
    parser.add_argument("--no-cap3",dest="cap3",action="store_false",help="Do not run CAP3, which joins the output of the different velvet runs")
    parser.add_argument("--no-exonerate",dest="exonerate",action="store_false",help="Do not run the Exonerate step, which assembles full length CDS regions and proteins from each gene")
    parser.add_argument("--velvet-mode",dest="use_velvet",action="store_true",help="Backwards compatability for velvet mode. NOT RECOMMENDED, VELVET MAKES ERROR PRONE ASSEMBLIES!")
    parser.add_argument("--no-assemble",dest="assemble",action="store_false",help="Skip the SPAdes assembly stage.")
    
    parser.add_argument('-r',"--readfiles",nargs='+',help="One or more read files to start the pipeline. If exactly two are specified, will assume it is paired Illumina reads.",default=[])
    parser.add_argument('-b','--baitfile',help="FASTA file containing bait sequences for each gene. If there are multiple baits for a gene, the id must be of the form: >Taxon-geneName",default=None)
    
    parser.add_argument('--cpu',type=int,default=0,help="Limit the number of CPUs. Default is to use all cores available.")
    parser.add_argument('--evalue',type=float,default=1e-10,help="e-value threshold for blastx hits, default: %(default)s")
    parser.add_argument('--max_target_seqs',type=int,default=10,help='Max target seqs to save in blast search, default: %(default)s')
    parser.add_argument('--cov_cutoff',type=int,default=8,help="Coverage cutoff for velvetg. default: %(default)s")
    parser.add_argument('--ins_length',type=int,default=200,help="Insert length for velvetg. default: %(default)s")
    parser.add_argument("--kvals",nargs='+',help="Values of k for velvet assemblies. Velvet needs to be compiled to handle larger k-values! Default auto-dectection by SPAdes.",default=None)
    parser.add_argument("--thresh",type=int,help="Percent Identity Threshold for stitching together exonerate results. Default is 55, but increase this if you are worried about contaminant sequences.",default=65)
    parser.add_argument("--length_pct",help="Include an exonerate hit if it is at least as long as X percentage of the reference protein length. Default = 90%%",default=90,type=int)
    parser.add_argument("--depth_multiplier",help="Accept any full-length exonerate hit if it has a coverage depth X times the next best hit. Set to zero to not use depth. Default = 10",default=10,type=int)

    parser.add_argument('--prefix',help="Directory name for pipeline output, default is to use the FASTQ file name.",default=None)
    parser.add_argument("--timeout",help="Use GNU Parallel to kill long-running processes if they take longer than X percent of average.",default=0)
    
    parser.add_argument("--target",help="Use this target to align sequences for each gene. Other targets for that gene will be used only for read sorting. Can be a tab-delimited file (one gene per line) or a single sequence name",default=None)
    parser.add_argument("--unpaired",help="Include a single FASTQ file with unpaired reads along with the two paired read files",default=False)
    parser.add_argument("--exclude",help="Do not use any sequence with the specified string as a target sequence for exonerate. The sequence will be used for read sorting.",default=None)
    parser.set_defaults(check_depend=False,blast=True,distribute=True,velvet=False,cap3=False,assemble=True,use_velvet=False,exonerate=True)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    
    run_dir = os.path.realpath(os.path.split(sys.argv[0])[0])
    print(("HybPiper was called with these arguments:\n{}\n".format(" ".join(sys.argv))))
    
    #Check dependencies
    if args.check_depend:
        if check_dependencies():
            other_scripts = ["distribute_reads_to_targets.py","distribute_targets.py","exonerate_hits.py"]
            for script in other_scripts:
                if os.path.isfile(os.path.join(run_dir,script)):
                    pass
                else:
                    print(("ERROR: Script {} not found! Please make sure it is in the same directory as this one!".format(script)))
                    return
            print("Everything looks good!")
            return
        else:
            print("ERROR: One or more dependencies not found!")
            return

    if args.baitfile:
        baitfile = os.path.abspath(args.baitfile)
    else:
        parser.print_help()
        return
    readfiles = [os.path.abspath(x) for x in args.readfiles]    
    if args.unpaired:
        unpaired_readfile = os.path.abspath(args.unpaired)
    else:
        unpaired_readfile = None
    if len(args.readfiles) < 1:
        print("ERROR: Please specify readfiles with -r")
        return
    if not args.baitfile:
        print("ERROR: Please specify a FASTA file containing target sequences.")
        return
    
    #Generate directory
    basedir, basename = make_basename(args.readfiles,prefix=args.prefix)
    os.chdir(os.path.join(basedir,basename))
    
    #BWA
    if args.bwa:
        if args.blast:
            args.blast=False
            bamfile = bwa(readfiles,baitfile,basename,cpu=args.cpu)
            if args.unpaired:
                unpaired_bamfile = bwa(unpaired_readfile,baitfile,basename,cpu=args.cpu,unpaired=True)
            if not bamfile:
                print("ERROR: Something went wrong with the BWA step, exiting!")
                return
        else:
            bamfile = basename + ".bam"
    
    #BLAST
    if args.blast:
        if args.unpaired:
            unpaired_blastxfile = blastx(unpaired_readfile,baitfile,args.evalue,basename,cpu=args.cpu,max_target_seqs=args.max_target_seqs,unpaired=True)
        blastx_outputfile = blastx(readfiles,baitfile,args.evalue,basename,cpu=args.cpu,max_target_seqs=args.max_target_seqs)
        if not blastx_outputfile:
            print("ERROR: Something is wrong with the Blastx step, exiting!")
            return
    else:
        blastx_outputfile = basename+".blastx"
    #Distribute
    
    if args.distribute:
        pre_existing_fastas = glob.glob("./*/*_interleaved.fasta") + glob.glob("./*/*_unpaired.fasta")
        for fn in pre_existing_fastas:
            os.remove(fn)
        if args.bwa:
            exitcode = distribute_bwa(bamfile,readfiles,baitfile,run_dir,args.target,unpaired_readfile,args.exclude)
        else:
            exitcode=    distribute(blastx_outputfile,readfiles,baitfile,run_dir,args.target,unpaired_readfile,args.exclude)
        if exitcode:
            sys.exit(1)
    if len(readfiles) == 2:
    
        genes = [x for x in os.listdir(".") if os.path.isfile(os.path.join(x,x+"_interleaved.fasta"))]
    else:
        genes = [x for x in os.listdir(".") if os.path.isfile(os.path.join(x,x+"_unpaired.fasta"))]
    #genes = ["gene008"]
    print("READFILES:\n")
    print(readfiles)
    print(genes)
    print(os.listdir('.'))
    if len(genes) == 0:
        print("ERROR: No genes with BLAST hits! Exiting!")
        return
    
    if args.use_velvet:
        #Velvet
        if args.velvet:
            exitcode = velvet(genes,cov_cutoff=args.cov_cutoff,ins_length=args.ins_length,kvals=args.kvals,cpu=args.cpu)
            if exitcode:
                return
            
        #CAP3
        if args.cap3:
            exitcode = cap3(genes,cpu=args.cpu)
            if exitcode:
                return

    if args.assemble:
        if len(readfiles) == 1:    
            spades_genelist = spades(genes,run_dir,cov_cutoff=args.cov_cutoff,cpu=args.cpu,kvals=args.kvals,paired=False,timeout=args.timeout)
        elif len(readfiles) == 2:
            if unpaired_readfile:
                spades_genelist = spades(genes,run_dir,cov_cutoff=args.cov_cutoff,cpu=args.cpu,kvals=args.kvals,timeout=args.timeout,unpaired=True)
            else:
                spades_genelist = spades(genes,run_dir,cov_cutoff=args.cov_cutoff,cpu=args.cpu,kvals=args.kvals,timeout=args.timeout)

        else:
            print("ERROR: Please specify either one (unpaired) or two (paired) read files! Exiting!")
            return
        if not spades_genelist:
            print("ERROR: No genes had assembled contigs! Exiting!")
            return
    
    #Exonerate hits
    if args.exonerate:
        genes = [x.rstrip() for x in open(exonerate_genefilename).readlines()]
        exitcode = exonerate(genes,basename,run_dir,cpu=args.cpu,thresh=args.thresh,length_pct = args.length_pct,depth_multiplier=args.depth_multiplier,timeout=args.timeout)
        if exitcode:
            return
    
    sys.stderr.write("Generated sequences from {} genes!\n".format(len(open("genes_with_seqs.txt").readlines())))
    
    paralog_warnings = [x for x in os.listdir(".") if os.path.isfile(os.path.join(x,basename,"paralog_warning.txt"))]
    with open("genes_with_paralog_warnings.txt",'w') as pw:
        pw.write("\n".join(paralog_warnings))
    sys.stderr.write("WARNING: Potential paralogs detected for {} genes!".format(len(paralog_warnings)))

if __name__ == "__main__":main()
