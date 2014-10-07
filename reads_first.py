#!/usr/bin/env python

import argparse,os,sys,importlib,shutil,subprocess

helptext="""
This script is a wrapper around several scripts in the HybSeqPipeline.
It can check whether you have the appropriate dependencies available (see --check-depend).
It makes sure that the other scripts needed are in the same directory as this one.
Command line options are passed to the other executables.
Unless --prefix is set, output will be put within a directory named after your read files."""

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
	executables = 		["blastx",
						"exonerate",
						"velvetg",
						"velveth",
						"cap3",
						"parallel"]

	python_packages = 	["Bio"]
	
	everything_is_awesome = True
	for e in executables:
		e_loc = py_which(e)
		if e_loc:
			print "{} found at {}".format(e,e_loc)
		else:
			print "{} not found in your $PATH!".format(e)
			everything_is_awesome = False
			
	for p in python_packages:
		try:
			i = importlib.import_module(p)
			print "Package {} successfully loaded!".format(p)
		except ImportError:
			print "Package {} not found!".format(p)
			everything_is_awesome=False
	return everything_is_awesome


def blastx(readfiles,baitfile,evalue,basename,cpu=None,max_target_seqs=10):
	if os.path.isfile(baitfile):
		if os.path.isfile(os.path.split(baitfile)[0]+'.psq'):
			db_file = baitfile
		else:
			print "Making protein blastdb in current directory."
			if os.path.split(baitfile)[0]:
				shutil.copy(baitfile,'.')
			db_file = os.path.split(baitfile)[1]
			makeblastdb_cmd = "makeblastdb -dbtype prot -in {}".format(db_file)
			print makeblastdb_cmd
			exitcode = subprocess.call(makeblastdb_cmd,shell=True)
			if exitcode:
				return None
	else:
		print "Cannot find baitfile at: {}".format(baitfile)
		return None

	#Remove previous blast results if they exist (because we will be appending)
	if os.path.isfile(basename+".blastx"):
		os.remove(basename+".blastx")
	
	for read_file in readfiles:
	
		
	
		#Piping commands for Fastq -> FASTA	
		# Curly braces must be doubled within a formatted string.
		pipe_cmd = "cat {} |  awk '{{if(NR % 4 == 1 || NR % 4 == 2) {{sub(/@/, \">\"); print; }} }}'".format(read_file)
	
		blastx_command = "blastx -db {} -query - -evalue {} -outfmt 6 -max_target_seqs {}".format(db_file,evalue,max_target_seqs)
		if cpu:
			full_command = "time {} | parallel -j {} --block 200K --recstart '>' --pipe {} >> {}.blastx ".format(pipe_cmd,cpu,blastx_command,basename)
		else:
			full_command = "time {} | parallel --block 200K --recstart '>' --pipe {} >> {}.blastx ".format(pipe_cmd,blastx_command,basename)
		print full_command
		exitcode = subprocess.call(full_command,shell=True)
		if exitcode:
			#Concatenate the two blastfiles.
		
			return None
			
	return basename + '.blastx'

def distribute(blastx_outputfile,readfiles,baitfile,run_dir):
	#NEED TO ADD SOMETHING ABOUT DIRECTORIES HERE.
	#print run_dir
	read_cmd = "time python {} {} {}".format(os.path.join(run_dir,"distribute_reads_to_targets.py"),blastx_outputfile," ".join(readfiles))
	exitcode = subprocess.call(read_cmd,shell=True)
	if exitcode:
		print "ERROR: Something went wrong with distributing reads to gene directories."
		return exitcode
	target_cmd = "time python {} {} --blastx {}".format(os.path.join(run_dir,"distribute_targets.py"),baitfile,blastx_outputfile)
	exitcode = subprocess.call(target_cmd,shell=True)
	if exitcode:
		print "ERROR: Something went wrong distributing targets to gene directories."
		return exitcode
	return None

def make_basename(readfiles,prefix=None):
	"""Unless prefix is set, generate a directory based off the readfiles, using everything up to the first underscore."""
	basename = os.path.split(readfiles[0])[1].split('_')[0]
	if not os.path.exists(basename):
		os.makedirs(basename)
	return basename

def velvet(genes,cov_cutoff=5,ins_length=200,kvals = ["21","31","41","51","61"],cpu=None,paired=True):
	"""Use parallel to run velveth and velvetg on a set of k values on every gene with blastx hits from the previous steps."""
	if os.path.isfile('velveth.log'):
		os.remove('velveth.log')
	if os.path.isfile('velvetg.log'):
		os.remove('velvetg.log')
	if paired:
		if cpu:
			velveth_cmd = "time parallel  -j {} --eta velveth {{1}}/velvet{{2}} {{2}} -shortPaired {{1}}/{{1}}_interleaved.fasta '>>' velveth.log ::: {} ::: {}".format(cpu," ".join(genes)," ".join(kvals))
		else:
			velveth_cmd = "time parallel  --eta velveth {{1}}/velvet{{2}} {{2}} -shortPaired {{1}}/{{1}}_interleaved.fasta '>>' velveth.log ::: {} ::: {}".format(" ".join(genes)," ".join(kvals))
	else:
		if cpu:
			velveth_cmd = "time parallel  -j {} --eta velveth {{1}}/velvet{{2}} {{2}} -short {{1}}/{{1}}_interleaved.fasta '>>' velveth.log ::: {} ::: {}".format(cpu," ".join(genes)," ".join(kvals))
		else:
			velveth_cmd = "time parallel  --eta velveth {{1}}/velvet{{2}} {{2}} -short {{1}}/{{1}}_interleaved.fasta '>>' velveth.log ::: {} ::: {}".format(" ".join(genes)," ".join(kvals))
	

	print os.getcwd()
	print os.listdir(".")
	
	
	print "Running velveth on {} genes".format(len(genes))
	print velveth_cmd
	exitcode = subprocess.call(velveth_cmd,shell=True)
	if exitcode:
		print "ERROR: Something went wrong with velveth!"
		return exitcode
	if cpu:
		velvetg_cmd = "time parallel -j {} --eta velvetg {{1}}/velvet{{2}} -ins_length {} -cov_cutoff {} '>>' velvetg.log ::: {} ::: {}".format(cpu,ins_length,cov_cutoff," ".join(genes)," ".join(kvals))
	else:
		velvetg_cmd = "time parallel --eta velvetg {{1}}/velvet{{2}} -ins_length {} -cov_cutoff {} '>>' velvetg.log ::: {} ::: {}".format(ins_length,cov_cutoff," ".join(genes)," ".join(kvals))
	print "Running velvetg on {} genes".format(len(genes))
	print velvetg_cmd
	exitcode = subprocess.call(velvetg_cmd,shell=True)
	if exitcode:
		print "ERROR: Something went wrong with velvetg!"
		return exitcode
	return None
	
def cap3(genes,cpu=None):
	print "Concatenating velvet output"
	if cpu:
		cat_cmd = "time parallel -j {} cat {{1}}/*/contigs.fa '>' {{1}}/velvet_contigs.fa ::: {}".format(cpu," ".join(genes))
	else:
		cat_cmd = "time parallel cat {{1}}/*/contigs.fa '>' {{1}}/velvet_contigs.fa ::: {}".format(" ".join(genes))
	print cat_cmd
	exitcode = subprocess.call(cat_cmd,shell=True)
	if exitcode:
		print "ERROR: Something went wrong while concatenating the velvet output files."
		return exitcode
	#Check that the velvet output actually has data in it.
	genes = [x for x in genes if os.path.getsize(os.path.join(x,'velvet_contigs.fa')) > 0]
	
	if len(genes) == 0:
		print "No genes left! Exiting!"
		return 1
		
	print "Running Cap3"
	if cpu:
		cap3_cmd = "time parallel -j {} --eta cap3 {{1}}/velvet_contigs.fa -o 20 -p 99 '>' {{1}}/{{1}}_cap3.log ::: {}".format(cpu," ".join(genes))
	else:
		cap3_cmd = "time parallel --eta cap3 {{1}}/velvet_contigs.fa -o 20 -p 99 '>' {{1}}/{{1}}_cap3.log ::: {}".format(" ".join(genes))
	print cap3_cmd
	exitcode = subprocess.call(cap3_cmd,shell=True)
	if exitcode:
		print "ERROR: Something went wrong with CAP3!"
		return exitcode

	print "Joining CAP3 contigs and singletons"
	join_cmd = "time parallel cat {{1}}/velvet_contigs.fa.cap.contigs {{1}}/velvet_contigs.fa.cap.singlets '>' {{1}}/{{1}}_cap3ed.fa ::: {}".format(" ".join(genes))
	print join_cmd
	exitcode = subprocess.call(join_cmd,shell=True)
	if exitcode:
		print "ERROR: Something went wrong joining the CAP3 output files!"
		return exitcode
		
	return None	

def exonerate(genes,basename,run_dir,replace=True,cpu=None):
	#Check that each gene in genes actually has CAP3 output
	#cap3_sizes = [os.stat(os.path.join(x,x+"_cap3ed.fa")).st_size for x in genes]
	#print cap3_sizes
	if replace:
		for g in genes:
			if os.path.isdir(os.path.join(g,basename)):
				shutil.rmtree(os.path.join(g,basename))
	genes = [x for x in genes if os.stat(os.path.join(x,x+"_cap3ed.fa")).st_size > 0]
	if len(genes) == 0:
		print "ERROR: No genes recovered for {}!".format(basename)
		return 1
	
	print "Running Exonerate to generate sequences for {} genes".format(len(genes))
	if cpu:
		exonerate_cmd = "time parallel -j {} python {} {{}}/{{}}_baits.fasta {{}}/{{}}_cap3ed.fa --prefix {{}}/{} ::: {}".format(cpu,os.path.join(run_dir,"exonerate_hits.py"),basename," ".join(genes))
	else:
		exonerate_cmd = "time parallel python {} {{}}/{{}}_baits.fasta {{}}/{{}}_cap3ed.fa --prefix {{}}/{} ::: {}".format(os.path.join(run_dir,"exonerate_hits.py"),basename," ".join(genes))
	print exonerate_cmd
	exitcode = subprocess.call(exonerate_cmd,shell=True)
	if exitcode:
		print "ERROR: Something went wrong with Exonerate!"
		return exitcode
	return
			
def main():
	parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("--check-depend",dest='check_depend',help="Check for dependencies (executables and Python packages) and exit. May not work at all on Windows.",action='store_true')
	parser.add_argument("--no-blast",dest="blast",action="store_false",help="Do not run the blast step. Downstream steps will still depend on the *_all.blastx file. \nUseful for re-runnning assembly/exonerate steps with different options.")
	parser.add_argument("--no-distribute",dest="distribute",action="store_false",help="Do not distribute the reads and bait sequences to sub-directories.")
	parser.add_argument("--no-velvet",dest="velvet",action="store_false",help="Do not run the velvet stages (velveth and velvetg)")
	parser.add_argument("--no-cap3",dest="cap3",action="store_false",help="Do not run CAP3, which joins the output of the different velvet runs")
	parser.add_argument("--no-exonerate",dest="exonerate",action="store_false",help="Do not run the Exonerate step, which assembles full length CDS regions and proteins from each gene")
	parser.add_argument('-r',"--readfiles",nargs='+',help="One or more read files to start the pipeline. If exactly two are specified, will assume it is paired Illumina reads.",default=[])
	parser.add_argument('-b','--baitfile',help="FASTA file containing bait sequences for each gene. If there are multiple baits for a gene, the id must be of the form: >Taxon-geneName",default=None)
	
	parser.add_argument('--cpu',type=int,default=0,help="Limit the number of CPUs. Default is to use all cores available.")
	parser.add_argument('--evalue',type=float,default=1e-9,help="e-value threshold for blastx hits, default: %(default)s")
	parser.add_argument('--max_target_seqs',type=int,default=10,help='Max target seqs to save in blast search, default: %(default)s')
	parser.add_argument('--cov_cutoff',type=int,default=4,help="Coverage cutoff for velvetg. default: %(default)s")
	parser.add_argument('--ins_length',type=int,default=200,help="Insert length for velvetg. default: %(default)s")
	parser.add_argument("--kvals",nargs='+',help="Values of k for velvet assemblies. Velvet needs to be compiled to handle larger k-values! Default is 21,31,41,51, and 61.",default=["21","31","41","51","61"])

	parser.add_argument('--prefix',help="Directory name for pipeline output, default is to use the FASTQ file name.",default=None)
	
	parser.set_defaults(check_depend=False,blast=True,distribute=True,velvet=True,cap3=True,exonerate=True)
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)
	args = parser.parse_args()
	
	run_dir = os.path.split(sys.argv[0])[0]

	#Check dependencies
	if args.check_depend:
		if check_dependencies():
			other_scripts = ["distribute_reads_to_targets.py","distribute_targets.py","exonerate_hits.py"]
			for script in other_scripts:
				if os.path.isfile(os.path.join(run_dir,script)):
					pass
				else:
					print "ERROR: Script {} not found! Please make sure it is in the same directory as this one!".format(script)
					return
			print "Everything looks good!"
		else:
			print "ERROR: One or more dependencies not found!"
			return

	if args.baitfile:
		baitfile = os.path.abspath(args.baitfile)
	else:
		parser.print_help()
		return
	readfiles = [os.path.abspath(x) for x in args.readfiles]	

	if len(args.readfiles) < 1:
		print "ERROR: Please specify readfiles with -r"
		return
	if not args.baitfile:
		print "ERROR: Please specify a FASTA file containing bait sequences."
		return
	
	#Generate directory
	basename = make_basename(args.readfiles,prefix=args.prefix)
	os.chdir(basename)
	
	#BLAST
	if args.blast:
		blastx_outputfile = blastx(readfiles,baitfile,args.evalue,basename,cpu=args.cpu,max_target_seqs=args.max_target_seqs)
		if not blastx_outputfile:
			print "ERROR: Something is wrong with the Blastx step, exiting!"
			return
	else:
		blastx_outputfile = basename+".blastx"
	#Distribute
	
	if args.distribute:
		exitcode=	distribute(blastx_outputfile,readfiles,baitfile,run_dir)
		if exitcode:
			sys.exit(1)
	genes = [x for x in os.listdir(".") if os.path.isfile(os.path.join(x,x+"_interleaved.fasta"))]
	
	if len(genes) == 0:
		print "ERROR: No genes with BLAST hits! Exiting!"
		return
	
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
	
	genes = [x for x in genes if os.path.getsize(os.path.join(x,'velvet_contigs.fa')) > 0]
	#Exonerate hits
	if args.exonerate:
		exitcode = exonerate(genes,basename,run_dir,cpu=args.cpu)
		if exitcode:
			return
	
	
	


if __name__ == "__main__":main()