#!/usr/bin/python3.6

####################################################
###                PRIMARY SCRIPT                ###
###   IMPORT PIPELINE SCRIPTS AS FUNCTIONS AND   ###
###   PARSE THROUGH AnnotaPipeline.config FILE   ###
####################################################

# USAGE: python3 AnnotaPipelineMain.py --seq protein.fasta

# --- IMPORT PACKAGES ----------------------------------------------------------

import os
import sys
import subprocess
import pathlib
import argparse
import configparser
import logging
import shutil
from shutil import which
from Bio import SeqIO

# --- SCRIPT LOCATION ----------------------------------------------------------

script_pwd = pathlib.Path(sys.argv[0]).absolute()

# AnnotaPipeline location
pipeline_pwd = script_pwd.parents[0]

# initial directory
home_dir_pwd = pathlib.Path.cwd()

# --- PARSER ARGUMENTS ---------------------------------------------------------

parser = argparse.ArgumentParser(
	add_help=False,
	description='''
	AnnotaPipeline
	Input sequence with [-s] (modify other parameters in AnnotaPipeline.config)
	Make sure all required parameters are given
	More instructions are in AnnotaPipeline.config

	If you already have protein file, give it through the flag -p, this way, Augustus prediction will be
	ignored. If you have a gff file for this protein file, give it throuth -gff flag to get complete annotations.

''',
	epilog=""" >>>> -s and -p are mutually exclusive arguments <<<<

And shall these hopeful words bring love inside your heart...""",
	formatter_class=argparse.RawTextHelpFormatter
)

# requiredNamed = parser.add_argument_group('required arguments')
optionalNamed = parser.add_argument_group('optional arguments')

# mandatory arguments
#   type (default): string
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument(
	'-s', '--seq', dest='seq',
	metavar='[genomic_file]',
	help='input genomic sequence file',
)
group.add_argument(
	'-p', '--prot', dest='protein',
	metavar='[protein_file]',
	help='input protein sequence file'
)

# optional arguments
#   no argument: uses default
#   type (default): string
optionalNamed.add_argument(
	'-c', '--config', dest='annotaconfig',
	metavar='[AnnotaPipeline.config]',
	default=pipeline_pwd / "AnnotaPipeline.config",
	help='configuration file for AnnotaPipeline'
)

optionalNamed.add_argument(
	'-gff', dest='gff',
	metavar='gff_file.gff',
	help='Gff file of Protein file'
)

# custom [--help] argument
optionalNamed.add_argument(
	'-h', '-help', '--help',
	action='help',
	default=argparse.SUPPRESS,
	help='What kind of sorcery is this?'
)

# arguments saved here
args = parser.parse_args()

# --- CREATE LogFile -----------------------------------------------------------

logger = logging.getLogger('AnnotaPipeline')


# ----------------------Redirect STDOUT and STDERR to logfile--------------


class StreamToLogger(object):
	"""
	Fake file-like stream object that redirects writes to a logger instance.
	"""

	def __init__(self, logger, log_level=logging.INFO):
		self.logger = logger
		self.log_level = log_level
		self.linebuf = ''

	def write(self, buf):
		temp_linebuf = self.linebuf + buf
		self.linebuf = ''
		for line in temp_linebuf.splitlines(True):
			if line[-1] == '\n':
				self.logger.log(self.log_level, line.rstrip())
			else:
				self.linebuf += line

	def flush(self):
		if self.linebuf != '':
			self.logger.log(self.log_level, self.linebuf.rstrip())
		self.linebuf = ''


logging.basicConfig(
	level=logging.DEBUG,
	format='%(asctime)s | %(name)s | %(levelname)s | %(message)s',
	datefmt='%d/%m/%Y %H:%M:%S',
	filename="AnnotaPipeline_LogFile.log",
	filemode='a'
)

stderr_logger = logging.getLogger('AnnotaPipeline')
sl = StreamToLogger(stderr_logger, logging.ERROR)
sys.stderr = sl

# --- PARSE ARGUMENTS FROM .config FILE ----------------------------------------

config_pwd = pathlib.Path(args.annotaconfig).absolute()

config = configparser.ConfigParser()  # call configparser
config.optionxform = str  # set default: get varibles as str
config.read(str(config_pwd))  # read configuration file
config.sections()  # get sections of each program


# to get variables type: config[_SECTION_].get('_variable_name_')

# --- FUNCTIONS ----------------------------------------------------------------


def check_parameters(sections):
	# Variables to check databases
	sp_verify = True
	nr_verify = True
	trembl_verify = True
	for section in sections:  # get box of variables
		for key in config[str(section)]:  # get variable for each box
			# check if it's not empty
			if args.protein is not None and key is config['AUGUSTUS']:
				pass
			else:
				if (len(config[str(section)].get(key))) < 1:
					# Check if any secondary database were given
					if str(section) == "EssentialParameters" and key == "specific_path_db":
						sp_verify = False
					elif str(section) == "EssentialParameters" and key == "nr_db_path":
						nr_verify = False
					elif str(section) == "EssentialParameters" and key == "trembl_db_path":
						trembl_verify = False
					else:
						logger.error(
							"Variable [" + key + "]"
							+ " from section [" + str(section) + "]"
							+ " is null"
						)
						logger.info("Exiting")
						exit(
							"Variable [" + key + "]"
							+ " from section [" + str(section) + "]"
							+ " is null"
						)
	# Exit and report error if there is more than one database or if there is no one
	if sum([sp_verify, nr_verify, trembl_verify]) == 0:
		logger.error("Error, there is no Secondary database, please review config file!")
		logger.info("Exiting")
		exit("Error, there is no Secondary database, please review config file!")
	if sum([sp_verify, nr_verify, trembl_verify]) == 2:
		logger.error("Error, there is two secondary databases, select one of them in config file")
		logger.info("Exiting")
		exit("Error, there is two secondary databases, select one of them in config file")


def fasta_fetcher(input_fasta, id_list, fetcher_output):
	wanted = set(id_list)
	records = (r for r in SeqIO.parse(input_fasta, "fasta") if r.id in wanted)
	count = SeqIO.write(records, fetcher_output, "fasta")
	if count < len(wanted):
		logger.info("IDs not found in input FASTA file")


def sequence_cleaner(fasta_file, min_length=0, por_n=100):
	output_file = open("Clear_" + fasta_file, "w+")
	# Using the Biopython fasta parse we can read our fasta input
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		# Take the current sequence
		sequence = str(seq_record.seq).upper()
		# Check if the current sequence is according to the user parameters
		if len(sequence) >= min_length and ((float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n):
			output_file.write(">" + str(seq_record.id) + "\n" + str(sequence) + "\n")
	output_file.close()


def check_file(file):
	if os.path.isfile(file) == 0:
		logging.error("File " + str(file) + " does not exist, please check earlier steps")
		logging.shutdown()
		sys.exit(1)
	elif os.path.getsize(file) == 0:
		logging.warning("File " + str(file) + " is empty, please check earlier steps")
		logging.warning("Trying to continue, but failures may occur")
		logger.info("Sometimes things don't work as we expect ...and that's ok. "
					"But, seriously, check your inputs and outputs and rerun")


def is_tool(name):
	if which(name) is None:
		logging.error("Program: " + str(name) + " must be avaliable in $PATH")
		logging.shutdown()
		sys.exit(1)


def augustus_run():
	# create AUGUSTUS directory

	augustus_dir = pathlib.Path(augustus_main['augustus_path'])
	augustus_bin = augustus_dir / "bin" / "augustus"
	augustus_config = augustus_dir / "config"

	logger.info("AUGUSTUS prediction has started")

	# AUGUSTUS: command line
	aug_config = str("--AUGUSTUS_CONFIG_PATH=" + str(augustus_config))
	aug_command = str(augustus_bin) + " " + str(aug_config)

	for variable in config['AUGUSTUS']:
		if variable != "augustus_path":
			# Check if it's a flag
			if str(interpro.get(variable)).lower() == "flag":
				aug_command += (
						" -" + str(variable))
			else:
				aug_command += (
						" --" + str(variable) + "=" + str(augustus_main.get(variable))
				)

	aug_command += (
			" " + str(seq_file)
			+ " > " + "AUGUSTUS_" + str(AnnotaBasename) + ".gff"
	)

	logger.info(str(aug_command))

	subprocess.getoutput(aug_command)

	# Check if expected file exists
	check_file(str("AUGUSTUS_" + str(AnnotaBasename) + ".gff"))

	logger.info("AUGUSTUS prediction is finished")

	# Parsing AUGUSTUS files
	logger.info("AUGUSTUS parsing has started")

	augustus_script = augustus_dir / "scripts" / "getAnnoFasta.pl"

	aug_file = "AUGUSTUS_" + str(AnnotaBasename)

	subprocess.run([
		"perl",
		str(augustus_script),
		str(aug_file + ".gff"),
		str("--seqfile=" + str(seq_file))
	]
	)

	logger.info("AUGUSTUS parsing is finished")

	# Declare as global, so you can modify out of fuction
	global aug_parsing
	aug_parsing = aug_file + ".aa"


def gfftofasta():
	subprocess.run([
		"python3",
		str(pipeline_pwd / "gfftofasta_parser.py"),
		"-gff",
		str(gff_file),
		"-annot",
		str("All_Annotated_Products.txt"),
		"-b",
		str(AnnotaBasename),
		"-faf",
		str(augustus_folder / str("Clear_" + aug_parsing)),
		"-org",
		str('"%s"' % str(AnnotaPipeline.get('organism')))
	]
	)


def run_fasta_to_GFF():
	subprocess.run([
		"python3",
		str(pipeline_pwd / "fasta_to_GFF.py"),
		"-gff",
		str(gff_file),
		"-all",
		str("All_Annotated_Products.txt"),
		"-b",
		str(AnnotaBasename)
	]
	)


# --- CHECK EACH BOX OF VARIABLES ----------------------------------------------
is_tool("blastp")
is_tool("perl")
is_tool("rpsblast")
is_tool("sed")

sections_config = config.sections()
check_parameters(sections_config)

# --- PREPARING SOME VARIABLES -------------------------------------------------
# \\ Check if user pass protein and gff file -> if it is, redirect variables
if args.seq is not None:
	seq_file = pathlib.Path(args.seq).absolute()
if args.protein is not None:
	prot_path = pathlib.Path(args.protein).absolute()
if args.gff is not None:
	gff_path = pathlib.Path(args.gff).absolute()

# Parameters from config file, each line is one script/software configuration
AnnotaPipeline = config['EssentialParameters']
AnnotaBasename = AnnotaPipeline['basename']
keyword_list = config['EssentialParameters']['key_words']
augustus_main = config['AUGUSTUS']
seq_cleaner = config['SequenceCleaner']
cdhit = config['CD-Hit']
interpro = config['INTERPROSCAN']
hmmscan = config['HMMSCAN']
blast = config['BLAST']
rpsblast = config['RPSBLAST']

# AnnotaPipeline main directory
home_dir = "AnnotaPipeline_" + str(AnnotaBasename)

# FileExistsError exception will be ignored
pathlib.Path(home_dir).mkdir(exist_ok=True)
os.chdir(home_dir)
annota_pwd = pathlib.Path(home_dir_pwd / home_dir)

# --- AUGUSTUS -----------------------------------------------------------------

augustus_folder = pathlib.Path(annota_pwd / str("1_GenePrediction_" + AnnotaBasename))
pathlib.Path(augustus_folder).mkdir(exist_ok=True)
os.chdir(augustus_folder)

# ===========================================================================
# Run augustus or start with protein file?
if args.protein is None:
	augustus_run()
else:
	# Copy protein file to AUGUSTUS path and padronize variable to run Annotapipeline after augustus
	shutil.copy2(prot_path, augustus_folder)
	aug_parsing = args.protein

# ============================================================================
# --- SEQUENCE CLEANER -------------------------------------------------------

logger.info("SEQUENCE CLEANER has started")

# Clean only with min_size
sequence_cleaner(str(aug_parsing), int(seq_cleaner.get('minsize_seq')))

# Check if expected file exists
check_file(str("Clear_" + aug_parsing))

logger.info("SEQUENCE CLEANER is finished. Please check Clear_" + aug_parsing)

os.chdir(annota_pwd)

# --- BLAST -----------------------------------------------------------------

# Create BLAST directory
blast_folder = pathlib.Path(annota_pwd / str("2_SimilarityAnalysis_" + AnnotaBasename))
pathlib.Path(blast_folder).mkdir(exist_ok=True)
os.chdir(blast_folder)

# BLAST: command line
logger.info("BLAST execution and parsing has started")

# Select secondary database from config file
if (len(AnnotaPipeline.get("specific_path_db"))) > 1:
	spdb_path = AnnotaPipeline.get("specific_path_db")
	flag_spdb = "-spdb"
elif (len(AnnotaPipeline.get("nr_db_path"))) > 1:
	spdb_path = AnnotaPipeline.get("nr_db_path")
	flag_spdb = "-nr"
elif (len(AnnotaPipeline.get("trembl_db_path"))) > 1:
	spdb_path = AnnotaPipeline.get("trembl_db_path")
	flag_spdb = "-trbl"

subprocess.run([
	"python3",
	str(str(pipeline_pwd / "blastp_parser.py")),
	"-s",
	str(augustus_folder / str("Clear_" + aug_parsing)),
	"-sp",
	str(AnnotaPipeline.get('swissprot_path_db')),
	"-basename",
	str(AnnotaBasename),
	str(flag_spdb),  # Flag for databse
	str(spdb_path),  # path to database
	"-id",
	str(blast.get('identity')),
	"-pos",
	str(blast.get('positivity')),
	"-cov",
	str(blast.get('coverage')),
	"-kw",
	str(keyword_list),
	"-t",
	str(AnnotaPipeline.get('threads'))
]
)

logger.info("BLAST execution and parsing is finished")

# Check if expected file exists
check_file(str(AnnotaBasename + "_no_hit_products.txt"))
check_file(str(AnnotaBasename + "_hypothetical_products.txt"))
check_file(str(AnnotaBasename + "_annotated_products.txt"))

os.chdir(annota_pwd)

# --- INTERPROSCAN -------------------------------------------------------

# Create INTERPROSCAN directory
interpro_folder = pathlib.Path(annota_pwd / str("3_FunctionalAnnotation_" + AnnotaBasename))
pathlib.Path(interpro_folder).mkdir(exist_ok=True)
os.chdir(interpro_folder)

# Preparing file that will be used by InteProScan
logger.info("Preparing file for InterProScan Hypothetical Proteins execution")

hypothetical_id = str(blast_folder / str(AnnotaBasename + "_hypothetical_products.txt"))
no_hit_id = str(blast_folder / str(AnnotaBasename + "_no_hit_products.txt"))
data1 = [line.strip() for line in open(hypothetical_id, "r")]
data2 = [line.strip() for line in open(no_hit_id, "r")]

fasta_fetcher(str(augustus_folder / str("Clear_" + aug_parsing)),
			  (data1 + data2),
			  "Hypothetical_Products.fasta")

# Check if expected file exists
check_file("Hypothetical_Products.fasta")

logger.info("InterProScan Hypothetical Proteins file preparation complete")

logger.info("Running INTERPROSCAN with Hypothetical Proteins")
# INTERPROSCAN: commandline

# General
interpro_command_line = str(AnnotaPipeline.get('interpro_exe')) + " -i Hypothetical_Products.fasta -o " \
						+ str(AnnotaBasename + "_interproscan_hypothetical_output.gff3") + \
						" -f GFF3 -t p -goterms -iprlookup"

# Optionals
for variable in config['INTERPROSCAN']:
	if str(interpro.get(variable)).lower() == "flag":
		interpro_command_line += (
				" -" + str(variable))
	else:
		interpro_command_line += (
				" -" + str(variable) + " " + str(interpro.get(variable))
		)

logger.info(str(interpro_command_line))

subprocess.getoutput(interpro_command_line)

# INTERPROSCAN parser (info_parser.py) can run without this result, but must be a valid file.
if os.path.isfile(str(AnnotaBasename + "_interproscan_hypothetical_output.gff3")) == 0:
	# Generate valid file
	subprocess.run(["touch", str(AnnotaBasename + "_interproscan_hypothetical_output.gff3")])
	logger.info("Interproscan analysis return no results, moving on without this results.")
	logger.warning("Check if your sequences have special characters (like *), remove it and rerun")

logger.info("INTERPROSCAN finished for hypothetical proteins")
logger.info("Preparing file for InterProScan Annotated Proteins execution")

annotated_file = str(blast_folder / str(AnnotaBasename + "_annotated_products.txt"))
os.system("cat " + annotated_file + " | cut -f 1 > Temp_annotated_products.txt")  # Using only IDs from the file
annotated_id = [line.strip() for line in open("Temp_annotated_products.txt", "r")]

fasta_fetcher(str(augustus_folder / str("Clear_" + aug_parsing)),
			  annotated_id,
			  "Annotated_Products.fasta")

# Check if expected file exists
check_file("Annotated_Products.fasta")

logger.info("InterProScan Annotated Proteins file preparation complete")
os.remove("Temp_annotated_products.txt")

logger.info("Running INTERPROSCAN with Annotated Proteins")
# INTERPROSCAN: commandline

# General
interpro_command_line = str(AnnotaPipeline.get('interpro_exe')) + " -i Annotated_Products.fasta -o " \
						+ str(AnnotaBasename + "_interproscan_annotated_output.gff3") + \
						" -f GFF3 -t p -goterms -iprlookup"

# Optionals
for variable in config['INTERPROSCAN']:
	if str(interpro.get(variable)).lower() == "flag":
		interpro_command_line += (
				" -" + str(variable))
	else:
		interpro_command_line += (
				" -" + str(variable) + " " + str(interpro.get(variable))
		)

logger.info(str(interpro_command_line))

subprocess.getoutput(interpro_command_line)

# INTERPROSCAN parser (info_parser.py) can run without this result, but must be a valid file.
if os.path.isfile(str(AnnotaBasename + "_interproscan_annotated_output.gff3")) == 0:
	# Generate valid file
	subprocess.run(["touch", str(AnnotaBasename + "_interproscan_annotated_output.gff3")])
	logger.info("Interproscan analysis return no results, moving on without this results.")
	logger.warning("Check if your sequences have special characters (like *), remove it and rerun")

logger.info("INTERPROSCAN finished for annotated proteins")

# --------------------------------------------------------------------

logger.info("Running HMMSCAN with Hypothetical Proteins")

# General
hmmscan_command_line = str(AnnotaPipeline.get('hmm_exe')) + " --cpu " + str(
	AnnotaPipeline.get("threads")) + " --tblout " \
					   + str(AnnotaBasename + "_hmmscan_output.txt") + \
					   " --noali"

# Optionals
for variable in config['HMMSCAN']:
	if str(hmmscan.get(variable)).lower() == "flag":
		hmmscan_command_line += (
				" --" + str(variable))
	else:
		# These specific arguments are passed through '-'
		if any(arg == str(variable) for arg in ("E", "Z", "T")):
			hmmscan_command_line += (
					" -" + str(variable) + " " + str(hmmscan.get(variable))
			)
		# Everyone else are passed through '--'
		else:
			hmmscan_command_line += (
					" --" + str(variable) + " " + str(hmmscan.get(variable))
			)

hmmscan_command_line += " " + str(AnnotaPipeline.get('pfam')) + \
						" Hypothetical_Products.fasta 1>/dev/null 2>hmmscan.err"

logger.info(str(hmmscan_command_line))

subprocess.getoutput(hmmscan_command_line)

# Check if expected file exists
check_file(str(AnnotaBasename + "_hmmscan_output.txt"))

logger.info("HMMSCAN finished")

# -------------------------------------------------------------------

logger.info("Running RPSBLAST with Hypothetical Proteins")

# General
rpsblast_command_line = "rpsblast -query Hypothetical_Products.fasta -out " \
						+ str(AnnotaBasename + "_rpsblast_output.outfmt6") + \
						" -db " + str(AnnotaPipeline.get('cdd')) + \
						" -outfmt \"6 qseqid sseqid sacc bitscore evalue ppos pident qcovs stitle\"" + \
						" -num_threads " + str(AnnotaPipeline.get("threads"))
# Optionals
for variable in config['RPSBLAST']:
	if str(rpsblast.get(variable)).lower() == "flag":
		rpsblast_command_line += (
				" -" + str(variable))
	else:
		rpsblast_command_line += (
				" -" + str(variable) + " " + str(rpsblast.get(variable))
		)

logger.info(str(rpsblast_command_line))

subprocess.getoutput(rpsblast_command_line)

# Check if expected file exists
check_file(str(AnnotaBasename + "_rpsblast_output.outfmt6"))

logger.info("RPSBLAST finished")

# -----------------------------------------------------------------------

logger.info("Parsing information from INTERPRO, HMMSCAN and RPSBLAST")

subprocess.run([
	"python3",
	str(pipeline_pwd / "funcannotation_parser.py"),
	"-interpro",
	str(AnnotaBasename + "_interproscan_hypothetical_output.gff3"),
	"-hmm",
	str(AnnotaBasename + "_hmmscan_output.txt"),
	"-rpsblast",
	str(AnnotaBasename + "_rpsblast_output.outfmt6"),
	"-basename",
	str(AnnotaBasename)
]
)

# Cleaning the house
try:
	os.remove("hmmscan.err")
	os.rmdir("temp/")
except:
	pass

logger.info("INTERPROSCAN, HMMSCAN and RPSBLAST execution and parsing is finished")

os.chdir(annota_pwd)

# -----------------------------------------------------------------------

logger.info("Preparing file for protein annotation")

subprocess.run([
	"python3",
	str(pipeline_pwd / "info_parser.py"),
	"-ipr1",
	str(interpro_folder / str(AnnotaBasename + "_interproscan_annotated_output.gff3")),
	"-ipr2",
	str(interpro_folder / str(AnnotaBasename + "_interproscan_hypothetical_output.gff3")),
	"-a",
	str(blast_folder / str(AnnotaBasename + "_annotated_products.txt")),
	"-hy",
	str(blast_folder / str(AnnotaBasename + "_hypothetical_products.txt")),
	"-nh",
	str(blast_folder / str(AnnotaBasename + "_no_hit_products.txt"))
]
)

# Check if expected file exists
check_file("All_annotation_products.txt")

try:
	# Sort annotations
	os.system("sort -V All_annotation_products.txt -o All_Annotated_Products.txt")
	os.remove("All_annotation_products.txt")
except:
	pass

logger.info("All_Annotated_Products.txt file is complete")

logger.info("Generating fasta file from All_Annotated_Products.txt")

# ------------  Defining what file will be used ---------------------------
if args.gff is not None and args.protein is not None:  # User gave protein file and gff file
	# Run parser to generate fasta_file
	gff_file = str(gff_path)
	gfftofasta()
	logger.info("Running fasta_to_GFF to get gff file annotated")
	run_fasta_to_GFF()
	logger.info("GFF file is ready - Check" + str(AnnotaBasename) + "_Annotated_GFF.gff")
elif args.protein is not None and args.gff is None:  # User gave only protein file
	logger.info("GFF file wasn't given, fasta file will have only annotations")
	logger.info("Running fasta_simple.py")
	subprocess.run([
		"python3",
		str(pipeline_pwd / "fasta_simple.py"),
		"-annot",
		str("All_Annotated_Products.txt"),
		"-b",
		str(AnnotaBasename),
		"-faf",
		str(augustus_folder / str("Clear_" + aug_parsing)),
		"-org",
		str('"%s"' % str(AnnotaPipeline.get('organism')))
	]
	)
	logger.info("GFF file wasn't given, skipping script fasta_to_GFF")
else:  # User selected run Augustus
	gff_file = augustus_folder / str("AUGUSTUS_" + str(AnnotaBasename) + ".gff")
	gfftofasta()
	logger.info("Running fasta_to_GFF to get gff file annotated")
	run_fasta_to_GFF()
	logger.info("GFF file is ready")


# Remove double quotes from organism
try:
	subprocess.run([
		"sed",
		"-i",
		's/\"//g',
		str("Annota_" + str(AnnotaBasename) + ".fasta")
	]
	)
except:
	pass

logger.info("Annota annotated the annotations on the annoted file.")

# -----------------------------------------------------------------------
