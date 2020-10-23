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
import datetime
import argparse
import configparser
import logging
from Bio import SeqIO

# --- PARSER ARGUMENTS ---------------------------------------------------------

parser = argparse.ArgumentParser(
    add_help=False,
    description='''
    AnnotaPipeline

    Input sequence with [-s] (modify other parameters in AnnotaPipeline.config)
    Make sure all required parameters are given
    More instructions are in AnnotaPipeline.config
    ''',
    epilog="""And shall these hopeful words bring love inside your heart...""",
    formatter_class=argparse.RawTextHelpFormatter
)

requiredNamed = parser.add_argument_group('required arguments')
optionalNamed = parser.add_argument_group('optional arguments')

# mandatory arguments
#   type (default): string
requiredNamed.add_argument(
    '-s', '--seq', dest='seq',
    metavar='[protein_file]',
    help='input sequence file',
    required=True
)

optionalNamed.add_argument(
    '-gff', dest='gff',
    metavar='[protein_file]',
    default=0,
    help='input gff to change annotations in this file'
)

# optional arguments
#   no argument: uses default
#   type (default): string
optionalNamed.add_argument(
    '-c', '--config', dest='annotaconfig',
    metavar='[AnnotaPipeline.config]',
    default="AnnotaPipeline.config",
    help=('configuration file for AnnotaPipeline')
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


# to get variables type: config[_SECTION_].get('_variable_name_')

# --- FUNCTIONS ----------------------------------------------------------------


def check_parameters(sections):
    for section in sections:  # get box of variables
        for key in config[str(section)]:  # get variable for each box
            # check if it's not empty
            if (len(config[str(section)].get(key))) < 1:
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
        if len(sequence) >= min_length and ((float(sequence.count("N"))/float(len(sequence)))*100 <= por_n):
            output_file.write(">" + str(seq_record.id) + "\n" + str(sequence) + "\n")
    output_file.close()



config = configparser.ConfigParser()  # call configparser
config.optionxform = str  # set default: get varibles as str
config.read(args.annotaconfig)  # read configuration file
config.sections()  # get sections of each program


# to get variables type: config[_SECTION_].get('_variable_name_')



# --- CHECK EACH BOX OF VARIABLES ----------------------------------------------

sections_config = config.sections()
check_parameters(sections_config)

# --- PREPARING SOME VARIABLES -------------------------------------------------

AnnotaPipeline = config['EssentialParameters']
AnnotaBasename = AnnotaPipeline['basename']
keyword_list = config['EssentialParameters']['key_words']
augustus_main = config['AUGUSTUS']
seq_cleaner = config['SequenceCleaner']
## cdhit = config['CD-Hit']
blast = config['BLAST']
interpro = config['INTERPROSCAN']
hmmscan = config['HMMSCAN']
rpsblast = config['RPSBLAST']


# --- INTERPROSCAN -------------------------------------------------------


# Preparing file that will be used by InteProScan
logger.info("Preparing file for InterProScan Hypothetical Proteins execution")

hypothetical_id = str(AnnotaBasename + "_hypothetical_products.txt")
no_hit_id = str(AnnotaBasename + "_no_hit_products.txt")
data1 = [line.strip() for line in open(hypothetical_id, "r")]
data2 = [line.strip() for line in open(no_hit_id, "r")]

fasta_fetcher(str(args.seq),
              (data1+data2),
              "Hypothetical_Products.fasta")

logger.info("InterProScan Hypothetical Proteins file preparation complete")

logger.info("Running INTERPROSCAN with Hypothetical Proteins")
# INTERPROSCAN: commandline

# General
interpro_command_line = "interproscan.sh -i Hypothetical_Products.fasta -o " \
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


run_info = True

# Check if there are result from interproscan

if os.path.isfile(str(AnnotaBasename + "_interproscan_hypothetical_output.gff3")) == 0:
	os.system("touch " + str(AnnotaBasename + "_interproscan_hypothetical_output.gff3"))
	logger.warn("Interproscan give no results for hyphotetical proteins")
	run_info = False


logger.info(str(interpro_command_line))

subprocess.getoutput(interpro_command_line)

logger.info("Preparing file for InterProScan Annotated Proteins execution")

annotated_file = str(str(AnnotaBasename + "_annotated_products.txt"))
os.system("cat " + annotated_file + " | cut -f 1 > Temp_annotated_products.txt") ### Using only IDs from the file
annotated_id = [line.strip() for line in open("Temp_annotated_products.txt", "r")]

fasta_fetcher(str(args.seq),
              annotated_id,
              "Annotated_Products.fasta")

logger.info("InterProScan Annotated Proteins file preparation complete")
#os.remove("Temp_annotated_products.txt")

logger.info("Running INTERPROSCAN with Annotated Proteins")
# INTERPROSCAN: commandline

# General
interpro_command_line = "interproscan.sh -i Annotated_Products.fasta -o " \
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

# Check if there are result from interproscan
if os.path.isfile(str(AnnotaBasename + "_interproscan_hypothetical_output.gff3")) == 0:
        os.system("touch " + str(AnnotaBasename + "_interproscan_hypothetical_output.gff3"))
        logger.warn("Interproscan give no results for annotated proteins")
        run_info = False


logger.info(str(interpro_command_line))

subprocess.getoutput(interpro_command_line)

logger.info("INTERPROSCAN finished")

# --------------------------------------------------------------------

logger.info("Running HMMSCAN with Hypothetical Proteins")

# General
hmmscan_command_line = "hmmscan --cpu " + str(AnnotaPipeline.get("threads")) + " --tblout " \
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

logger.info("RPSBLAST finished")

# -----------------------------------------------------------------------

logger.info("Parsing information from INTERPRO, HMMSCAN and RPSBLAST")

subprocess.run([
    "python3",
    str("funcannotation_parser.py"),
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
if run_info == True:
	os.remove("Annotated_Products.fasta")
else:
	os.remove("Annotated_Products.fasta")
	os.remove("Hypothetical_Products.fasta")
	os.remove("hmmscan.err")
	os.rmdir("temp/")

logger.info("INTERPROSCAN, HMMSCAN and RPSBLAST execution and parsing is finished")


# -----------------------------------------------------------------------

if run_info == True:
	logger.info("Preparing file for protein annotation")

	subprocess.run([
    		"python3",
    		str("info_parser.py"),
    		"-ipr1",
    		str(AnnotaBasename + "_interproscan_annotated_output.gff3"),
    		"-ipr2",
    		str(AnnotaBasename + "_interproscan_hypothetical_output.gff3"),
    		"-a",
    		str(AnnotaBasename + "_annotated_products.txt"),
    		"-hy",
    		str(AnnotaBasename + "_hypothetical_products.txt"),
    		"-nh",
    		str(AnnotaBasename + "_no_hit_products.txt")
	]
	)

	logger.info("All_Annotated_Products.txt file is complete")

	logger.info("Generating fasta file from All_Annotated_Products.txt")

	if args.gff != 0:
		try:
    			subprocess.run([
    				"python3",
    				str("gfftofasta_parser.py"),
    				"-gff",
    				str(args.gff),
    				"-annot",
    				str("All_annotation_products.txt"),
    				"-b",
    				str(AnnotaBasename),
    				"-faf",
    				str(args.seq),
    				"-org",
    				str('"%s"' % str(AnnotaPipeline.get('organism')))
    			]
    			)
		except:
    			logger.info("Cannot reannote gff file --> skipping")
		pass

		os.system("sort -V All_annotation_products.txt -o All_Annotated_Products.txt")
		os.remove("All_annotation_products.txt")

logger.info("Annota annotated the annotations on the annoted file.")

# -----------------------------------------------------------------------
