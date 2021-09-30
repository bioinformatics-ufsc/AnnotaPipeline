#!/usr/bin/python3.6

####################################################
###                PRIMARY SCRIPT                ###
###   IMPORT PIPELINE SCRIPTS AS FUNCTIONS AND   ###
###   PARSE THROUGH AnnotaPipeline.config FILE   ###
####################################################

# USAGE: python3 AnnotaPipelineMain.py --seq protein.fasta

# --- IMPORT PACKAGES ----------------------------------------------------------

from Bio import SeqIO
from shutil import which
import pandas as pd
import argparse
import configparser
import logging
import os
import pathlib
import shutil
import subprocess
import sys
import re

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

   _           _                     _    __       _             _      
  | |         | |                   | |  /_/      (_)           | |     
  | |     __ _| |__   ___  _ __ __ _| |_ ___  _ __ _  ___     __| | ___ 
  | |    / _` | '_ \ / _ \| '__/ _` | __/ _ \| '__| |/ _ \   / _` |/ _ \
  | |___| (_| | |_) | (_) | | | (_| | || (_) | |  | | (_) | | (_| |  __/
  |______\__,_|_.__/ \___/|_|  \__,_|\__\___/|_|  |_|\___/   \__,_|\___|
                                                                        
                                                                        
   ____  _       _        __                       __  _   _            
  |  _ \(_)     (_)      / _|                     /_/ | | (_)           
  | |_) |_  ___  _ _ __ | |_ ___  _ __ _ __ ___   __ _| |_ _  ___ __ _  
  |  _ <| |/ _ \| | '_ \|  _/ _ \| '__| '_ ` _ \ / _` | __| |/ __/ _` | 
  | |_) | | (_) | | | | | || (_) | |  | | | | | | (_| | |_| | (_| (_| | 
  |____/|_|\___/|_|_| |_|_| \___/|_|  |_| |_| |_|\__,_|\__|_|\___\__,_| 
                                                                        
                                                                        
                      _    _ ______ _____  _____                        
                     | |  | |  ____/ ____|/ ____|                       
                     | |  | | |__ | (___ | |                            
                     | |  | |  __| \___ \| |                            
                     | |__| | |    ____) | |____                        
                      \____/|_|   |_____/ \_____|                       
                                                                        

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
    help='Gff file of protein file'
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

# CREATE CUSTOM LOG
# JONAS = 5 # Numeric error
# logging.addLevelName(5, "JONAS")
# logger.setLevel(JONAS)
# logger.log(JONAS,"teste")

# --- PARSE ARGUMENTS FROM .config FILE ----------------------------------------

config_pwd = pathlib.Path(args.annotaconfig).absolute()

config = configparser.ConfigParser()  # call configparser
config.optionxform = str  # set default: get varibles as str
config.read(str(config_pwd))  # read configuration file
config.sections()  # get sections of each program


# to get variables type: config[_SECTION_].get('_variable_name_')

# --- FUNCTIONS ----------------------------------------------------------------

def log_quit():
    logger.info("Exiting")
    logging.shutdown()
    sys.exit(1)


def is_tool(name):
    if which(name) is None:
        logging.error(f"Program: {str(name)} must be avaliable in $PATH")
        log_quit()


def kallisto_run(python_path, kallisto_exe, paired_end, method, basename, fasta):
    
    logger.info("KALLISTO index has started")
    kallisto_command_index = f"{kallisto_exe} index -i {basename}_kallisto_index.idx {fasta}"
    logger.debug(f"{kallisto_command_index}")
    subprocess.getoutput(kallisto_command_index)
    check_file(f"{basename}_kallisto_index.idx")

    # Standart command line for kallisto index
    # kallisto index -i transcripts.idx transcripts.fasta
    if paired_end == True:
        if len(kallisto.get("l")) != 0:
            l_flag = f"-l {kallisto.get('l')}"
            logger.info(f"KALLISTO will run with -l {kallisto.get('l')}")
        else:
            l_flag = ""
        if len(kallisto.get("s")) != 0:
            s_flag = f"-s {kallisto.get('s')}"
            logger.info(f"KALLISTO will run with -s {kallisto.get('s')}")
        else:
            s_flag = ""
       
        kallisto_command_quant = (
            f"{kallisto_exe} quant -i {basename}_kallisto_index.idx "
            f"{s_flag} {l_flag} "
            f"-o {basename}_kallisto_output -b {kallisto.get('bootstrap')} {kallisto.get('rnaseq-data')}"
        )
        logger.info(f"KALLISTO quant has started")
        logger.debug(kallisto_command_quant)
        subprocess.getoutput(kallisto_command_quant)
        # Standart command line for kallisto quant paired end
        # kallisto quant -i transcripts.idx -o output -b 100 reads_1.fastq reads_2.fastq
    else:
        kallisto_command_quant = (
            f"{kallisto_exe} quant -i {basename}_kallisto_index.idx "
            f"-l {kallisto.get('s')} -s {kallisto.get('l')} --single "
            f"-o kallisto_output -b {kallisto.get('bootstrap')} {kallisto.get('rnaseq-data')}"
        )
        logger.info(f"KALLISTO quant has started")
        logger.debug(kallisto_command_quant)
        subprocess.getoutput(kallisto_command_quant)
        # Standart command line for kallisto quant single end
        # kallisto quant -i transcripts.idx -o output -b 100 --single -l 180 -s 20 reads_1.fastq
    check_file(f"{basename}_kallisto_output/abundance.tsv")
    # Define method to parse
    if method == "median":
        kallisto_parser_flag = "-tpmmd"
    elif method == "mean":
        kallisto_parser_flag = "-tpmavg"
    else:
        kallisto_parser_flag = f'-tpmval {kallisto.get("value")}'
    # Run parser
    # inside 4_Transcript_Quantification_
    kallisto_parser_path =  str(pipeline_pwd / "kallisto_parser.py")
    kallisto_parser_command = (
        f"{python_path} {kallisto_parser_path} "
        f"-ktfile {basename}_kallisto_output/abundance.tsv " # kallisto output is always "abundance.tsv"
        f"-basename {basename} {kallisto_parser_flag}"
    )
    logger.info(f"KALLISTO parsing has started")
    subprocess.getoutput(kallisto_parser_command)
    check_file(f"{basename}_transcript_evidence.tsv")


def kallisto_check_parameters():
    # kallisto method indicate wich method to parse (mean, median, value) 
    # kallisto paired_end indicate rnaseq type (important to run kallisto)
    global kallisto_method, kallisto_paired_end
    # This variable will store method to parse kallisto output
    kallisto_method = None
    # if kallisto path were given, check other arguments, if not pass.
    if len(config[str('KALLISTO')].get("kallisto_path")) == 0:
        logger.info("Arguments for KALLISTO are empty. This step will be skipped.")
    else:
        kallisto_check = []
        for argument in ("rnaseq-data", "median", "mean", "value"):
            if argument == "rnaseq-data":
                if len(config[str('KALLISTO')].get("rnaseq-data").split()) > 2:
                    logger.error(
                        "[KALLISTO]: there are more arguments than required for rnaseq-data (KALLISTO). " \
                        "Pass one if your data is from single-end, " \
                        "and two files if your data is from paired-end!"
                    )
                    log_quit()
                ###### This box check how many files were given in rnaseq-data #####
                elif len(config['KALLISTO'].get("rnaseq-data").split()) == 2:
                    kallisto_paired_end = True
                    logger.info(f"KALLISTO will run with paired end data")
                    kallisto_check.append(argument)
                elif len(config['KALLISTO'].get("rnaseq-data").split()) == 1:
                    kallisto_paired_end = False
                    # Check required arguments for single end data
                    if len(config['KALLISTO'].get('l')) == 0:
                        logger.error("[KALLISTO]: mandatory argument for single-end data 'l' is empty")
                        log_quit()
                    elif len(config['KALLISTO'].get('s')) == 0:
                        logger.error("[KALLISTO]: mandatory argument for single-end data 's' is empty")
                        log_quit()
                    logger.info(f"KALLISTO will run with single-end data")
                    kallisto_check.append(argument)
                else:
                    logger.error(
                        "[KALLISTO]: check values for rnaseq-data (KALLISTO). " \
                        "Pass one if your data is from single-end, " \
                        "and two files if your data is from paired-end!"
                    )
                    log_quit()
                ####################################################################
            # Check if argument is Empty
            elif len(config[str('KALLISTO')].get(argument)) == 0:
                pass
            else:
                # If it's not empty, save argument
                kallisto_check.append(argument)
        # Check kallisto arguments
        # If all kallisto arguments are empty, then don't run this guy
        if len(kallisto_check) == 0:
            logger.info("Arguments for KALLISTO are empty. This step will be skipped.")
            kallisto_method = None
        # Rna-seq data is required for calisto
        if 'rnaseq-data' in kallisto_check:
            # if there is rna-seq data, check if method is correctly given
            if len(kallisto_check) > 2:
                logger.error("[KALLISTO]: there is more than one method selected to parse KALLISTO ouput. Please, review .config file.")
                log_quit()
            elif len(config[str('KALLISTO')].get("bootstrap")) == 0:
                logger.error("[KALLISTO]: bootstrap is empty, default value is 0. At least pass this value")
                log_quit()
            else:
                logger.info(f"KALLISTO will run with method: {kallisto_check[1]}")
                # Pass method, to use further
                kallisto_method = kallisto_check[1]                    


def percolator_check_parameters():
    if len(config['PERCOLATOR'].get('percolator_bash')) == 0:
        logger.error("[PERCOLATOR]: path to software is empty, but comet was setted")
        logger.warning("[PERCOLATOR]: leave comet fields empty or pass parameters for PERCOLATOR")
        log_quit()
    if len(config['PERCOLATOR'].get('qvalue')) == 0:
        logger.error("[PERCOLATOR]: qvalue cutoff is empty, check this parameter")
        log_quit()
    # check value for parser
    if not (0 <= float(config['PERCOLATOR'].get('qvalue')) <= 1):
        logger.error("[PERCOLATOR]: qvalue cutoff invalid. Must be float between [0-1]")
        log_quit()


def comet_check_parameters():
    if len(config['COMET'].get("comet_bash")) == 0:
        logger.info("Arguments for COMET are empty. This step will be skipped.")
    else:
        # ------------ check F and L -------------------
        last_check = False
        first_check = False
        if len(config['COMET'].get('first')) !=0:
            first_check = True
        if len(config['COMET'].get('last')) !=0:
            last_check = True
        # ----------------------------------------------
        for argument in ("params", "mass_files"):
            if len(config['COMET'].get(argument)) == 0:
                logger.error(f"[COMET]: Parameter [{argument}] from section [COMET] is null")
                log_quit()
        # ------------ check extension mass files ------
        if len(config['COMET'].get('mass_files_ext').split()) < 1:
            logger.error(f"[COMET]: Parameter [{config['COMET'].get('mass_files_ext')}] from section [COMET] is null")
            log_quit()
        elif len(config['COMET'].get('mass_files_ext').split()) > 1:
            logger.error(f"[COMET]: Parameter [{config['COMET'].get('mass_files_ext')}] from section [COMET] have more than one argument")
            log_quit()
        # ------------ Specific conditions -----------------------
        # XOR conditional >> check if one value is True and another is False
        if last_check ^ first_check:
            logger.error("[COMET]: both arguments from comet, first and last, must be given")
            logger.error("[COMET]: Leave both empty or give both")
            log_quit()
        elif last_check == first_check == True:
            logger.warning("[COMET]: values for 'fist' and 'last' will overwrite those in params file")
            global use_last_and_first
            use_last_and_first = True
        else:
            use_last_and_first = False
        # Only check percolator if comet parames are ok
        percolator_check_parameters()


def check_parameters(sections):
    # Variables to check databases
    # swissprot database
    sp_verify = True
    # nr database
    nr_verify = True
    # trembl database
    trembl_verify = True
    for section in sections:  # get box of variables
        # Check kallisto optional arguments
        if str(section) == "KALLISTO":
            kallisto_check_parameters()  
        elif str(section) == "COMET":
            comet_check_parameters()
        elif str(section) == "PERCOLATOR":
            # Percolator depends on COMET, thus, is checked with comet_check_parameters()
            pass
        else:  
            for key in config[str(section)]:  # get variable for each box
                # Arguments for agutustus don't need to be checked if protein file were given
                if args.protein is not None and key is config['AUGUSTUS']:
                    pass                       
                else:
                    if (len(config[str(section)].get(key))) < 1:
                        # Check if any secondary database were given
                        if str(section) == "EssentialParameters" and key == "specific_db_path":
                            sp_verify = False
                        elif str(section) == "EssentialParameters" and key == "nr_db_path":
                            nr_verify = False
                        elif str(section) == "EssentialParameters" and key == "trembl_db_path":
                            trembl_verify = False
                        # Any other parameter checked
                        else:
                            # Crash pipeline if some required variable is empty
                            logger.error(f"Variable [{key}] from section [{str(section)}] is null")
                            log_quit()
    # Exit and report error if there is more than one database or if there is no one
    if sp_verify == nr_verify == trembl_verify == False:
        logger.error("[DATABASE]: there is no secondary database. Please review config file!")
        log_quit()
    if sum([sp_verify, nr_verify, trembl_verify]) == 2:
        logger.error("[DATABASE]: there are two secondary databases. Select one of them in the config file.")
        log_quit()


def fasta_fetcher(input_fasta, id_list, fetcher_output):
    wanted = sorted(set(id_list))
    #records = (r for r in SeqIO.parse(input_fasta, "fasta") if r.id in wanted)
    records = (seq for seq in SeqIO.parse(input_fasta, "fasta") for r in wanted if r in seq.id)
    count = SeqIO.write(records, fetcher_output, "fasta")
    if count < len(wanted):
        logger.info("IDs not found in input FASTA file")


def annotate_codingseq(aa_fasta, codingseq_fasta, basename):
    # Linha de anotacao completa
    anno_all = [line.strip() for line in open(aa_fasta) if ">" in line]
    anno_all = [sub.replace(">", "") for sub in anno_all]
    # IDs simplificadas
    anno_ids = [name.split("|")[0].strip() for name in anno_all]

    corrrected_ids = dict(zip(anno_ids, anno_all))

    id_dict  = SeqIO.to_dict(SeqIO.parse(codingseq_fasta, "fasta"))

    corrected_fasta = str(f"AnnotaPipeline_{basename}_transcripts.fasta")

    logger.info(f"Generating AnnotaPipeline_{basename}_transcripts.fasta file from AnnotaPipeline_{basename}_proteins.fasta")

    with open(corrected_fasta, "w") as corrected:
        for key, record in id_dict.items():
            for id_key in corrrected_ids.keys():
                if id_key in key:
                    record.id = corrrected_ids.get(id_key)
                    record.description = ""
                    record.seq = record.seq.upper()
                    SeqIO.write(record, corrected, "fasta-2line")


def sequence_cleaner(fasta_file, min_length=0, por_n=100):
    output_file = open(f"Clear_{fasta_file}", "w+")
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        # Take the current sequence
        sequence = str(seq_record.seq).upper()
        # Check if the current sequence is according to the user parameters
        if len(sequence) >= min_length and ((float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n):
            output_file.write(f">{str(seq_record.id)}\n{str(sequence)}\n")
    output_file.close()


def check_file(file):
    logger = logging.getLogger('AnnotaPipeline')
    if os.path.isfile(file) == 0:
        logging.error(f"File {str(file)} does not exist, please check earlier steps")
        log_quit()
    elif os.path.getsize(file) == 0:
        logging.warning(f"File {str(file)} is empty, please check earlier steps")
        logging.warning("Trying to continue, but failures may occur")
        logger.info("Sometimes things don't work as we expect... and that's ok. "
                    "But, seriously, check your inputs and outputs and rerun.")


def augustus_run(basename):
    # create AUGUSTUS directory
    logger = logging.getLogger('AUGUSTUS')
    augustus_dir = pathlib.Path(augustus_main['augustus_path'])
    augustus_bin = augustus_dir / "bin" / "augustus"
    augustus_config = augustus_dir / "config"

    logger.info("AUGUSTUS prediction has started")

    # AUGUSTUS: command line
    aug_config = f"--AUGUSTUS_CONFIG_PATH={str(augustus_config)}"
    aug_command = f"{str(augustus_bin)} {str(aug_config)}"

    for variable in config['AUGUSTUS']:
        if variable != "augustus_path":
            # Check if it's a flag
            if str(interpro.get(variable)).lower() == "flag":
                aug_command += f" -{str(variable)}"
            else:
                aug_command += f" --{str(variable)}={str(augustus_main.get(variable))}"

    aug_command += f" {str(seq_file)} > AUGUSTUS_{str(basename)}.gff"

    logger.debug(str(aug_command))

    subprocess.getoutput(aug_command)

    # Check if expected file exists
    check_file(f"AUGUSTUS_{str(basename)}.gff")

    logger.info("AUGUSTUS prediction is finished")

    # Parsing AUGUSTUS files
    logger.info("AUGUSTUS parsing has started")

    augustus_script = augustus_dir / "scripts" / "getAnnoFasta.pl"

    aug_file = f"AUGUSTUS_{str(basename)}"

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


def gfftofasta(python_path):
    subprocess.run([
        python_path,
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


def run_fastatogff(python_path):
    subprocess.run([
        python_path,
        str(pipeline_pwd / "fastatogff.py"),
        "-gff",
        str(gff_file),
        "-all",
        str("All_Annotated_Products.txt"),
        "-b",
        str(AnnotaBasename)
        ]
    )


def add_features(feature, data, save):
    # Work on temp dataset
    parcial_feature = data.copy()
    # Count uniq entries and place new boolean column
    parcial_feature[f"Unique {feature}"] = parcial_feature.index.isin(parcial_feature.drop_duplicates(f"{feature}", keep =False).index)
    # Replace boolean values for 0 and 1 (easiest do count)
    parcial_feature.replace({False: 0, True: 1}, inplace=True)
    # Group total peptides/spectrum by proteinIndex
    parcial_feature_2 = parcial_feature.groupby(['ProteinID']).size().sort_values(ascending=False).reset_index(name=f'Total {feature}')
    # Extract rows with unique peptides/spectrum (assigned with 1)
    feature_process = parcial_feature.loc[parcial_feature[f"Unique {feature}"] == 1].drop(columns=[f"Peptide", "Spectrum"])
    # Count how many unique values are for each protein
    feature_process = feature_process.groupby(['ProteinID']).size().reset_index(name=f'Unique {feature}')
    if 'ProteinID' not in save.columns:
        save["ProteinID"] = parcial_feature_2["ProteinID"]
    # Add unique feature column
    save = save.set_index("ProteinID").join(feature_process.set_index("ProteinID")).reset_index()
    # Add total feature column
    save = save.set_index("ProteinID").join(parcial_feature_2.set_index("ProteinID")).reset_index()
    return save



def quantitative_proteomics(path, basename):
    # get all percolator parsed files
    parsed_files = pathlib.Path(path).glob('*_parsed.tsv')

    # Start Empty dataframe to store all _parsed files
    data = pd.DataFrame({'ProteinID':[], \
                        'Peptide':[], \
                        'Spectrum':[]}) 
    for file in parsed_files:
        df = pd.read_csv(file, sep='\t', header=0)
        data = data.append(df)

    # Empty dataframe to save
    total = pd.DataFrame({}) 
    total = add_features('Peptide', data, total)
    total = add_features('Spectrum', data, total)
    total = total.fillna(0).astype({"Unique Peptide": int, "Unique Spectrum": int}).sort_values(by='ProteinID', ascending=False)
    total.to_csv(f"{basename}_pre_total_Proteomics_Quantification.tsv", sep="\t", index=False)


# ------------------------------------------------------------------------------
# --- CHECK EACH BOX OF VARIABLES ----------------------------------------------
is_tool("blastp")
is_tool("perl")
is_tool("rpsblast")
# ------------------------------------------------------------------------------

logger.info("")

sections_config = config.sections()
logger = logging.getLogger('AnnotaPipeline')

logger.info("Checking parameters in AnnotaPipeline.config file")
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
python_exe = config['EssentialParameters']['python_exe']
AnnotaPipeline = config['EssentialParameters']
AnnotaBasename = AnnotaPipeline['basename']
keyword_list = config['EssentialParameters']['key_words']
augustus_main = config['AUGUSTUS']
seq_cleaner = config['SequenceCleaner']
interpro = config['INTERPROSCAN']
hmmscan = config['HMMSCAN']
blast = config['BLAST']
rpsblast = config['RPSBLAST']
kallisto = config['KALLISTO']
comet = config['COMET']
percolator = config['PERCOLATOR']

# AnnotaPipeline main directory
home_dir = f"AnnotaPipeline_{str(AnnotaBasename)}"

# FileExistsError exception will be ignored
pathlib.Path(home_dir).mkdir(exist_ok=True)
os.chdir(home_dir)
annota_pwd = pathlib.Path(home_dir_pwd / home_dir)

# --- AUGUSTUS -----------------------------------------------------------------

augustus_folder = pathlib.Path(annota_pwd / str("1_GenePrediction_" + AnnotaBasename))
pathlib.Path(augustus_folder).mkdir(exist_ok=True)
os.chdir(augustus_folder)

# ==============================================================================
# Run Augustus or start with protein file?
if args.protein is None:
    augustus_run(AnnotaBasename)
else:
    # Copy protein file to AUGUSTUS path and padronize variable to run Annotapipeline after augustus
    shutil.copy2(prot_path, augustus_folder)
    aug_parsing = args.protein

# ==============================================================================
# SEQUENCE CLEANER -------------------------------------------------------------
logger = logging.getLogger('AnnotaPipeline')
logger.info("Sequence Cleaner has started")

# Clean only with min_size
sequence_cleaner(str(aug_parsing), int(seq_cleaner.get('minsize_seq')))

# Check if expected file exists
check_file(f"Clear_{aug_parsing}")

logger.info(f"Sequence Cleaner is finished. Please check Clear_{aug_parsing}")

os.chdir(annota_pwd)

# BLAST ------------------------------------------------------------------------

# Create BLAST directory
blast_folder = pathlib.Path(annota_pwd / str("2_SimilarityAnalysis_" + AnnotaBasename))
pathlib.Path(blast_folder).mkdir(exist_ok=True)
os.chdir(blast_folder)

# BLAST: command line
logger.info("BLAST execution and parsing has started")

# Select secondary database from config file
if (len(AnnotaPipeline.get("specific_db_path"))) > 1:
    spdb_path = AnnotaPipeline.get("specific_db_path")
    flag_spdb = "-spdb"
elif (len(AnnotaPipeline.get("nr_db_path"))) > 1:
    spdb_path = AnnotaPipeline.get("nr_db_path")
    flag_spdb = "-nr"
elif (len(AnnotaPipeline.get("trembl_db_path"))) > 1:
    spdb_path = AnnotaPipeline.get("trembl_db_path")
    flag_spdb = "-trbl"

subprocess.run([
    str(python_exe),
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
check_file(f"{str(AnnotaBasename)}_no_hit_products.txt")
check_file(f"{str(AnnotaBasename)}_hypothetical_products.txt")
check_file(f"{str(AnnotaBasename)}_annotated_products.txt")

os.chdir(annota_pwd)

# INTERPROSCAN -----------------------------------------------------------------

# Create INTERPROSCAN directory
interpro_folder = pathlib.Path(annota_pwd / str("3_FunctionalAnnotation_" + AnnotaBasename))
pathlib.Path(interpro_folder).mkdir(exist_ok=True)
os.chdir(interpro_folder)

logger = logging.getLogger('INTERPROSCAN')
# Preparing file that will be used by InteProScan
logger.info("Preparing file for INTERPROSCAN Hypothetical Proteins execution")

hypothetical_id = str(blast_folder / str(AnnotaBasename + "_hypothetical_products.txt"))
no_hit_id = str(blast_folder / str(AnnotaBasename + "_no_hit_products.txt"))
hypothetical_id_strip = [line.strip() for line in open(hypothetical_id, "r")]
no_hit_id_strip = [line.strip() for line in open(no_hit_id, "r")]

fasta_fetcher(str(augustus_folder / str("Clear_" + aug_parsing)),
              (hypothetical_id_strip + no_hit_id_strip),
              "Hypothetical_Products.fasta")

# Check if expected file exists
check_file("Hypothetical_Products.fasta")

logger.info("INTERPROSCAN Hypothetical Proteins file preparation complete")

logger.info("Running INTERPROSCAN with Hypothetical Proteins")
# INTERPROSCAN: commandline

# General
interpro_command_line = (
    f"{str(AnnotaPipeline.get('interpro_exe'))} -i Hypothetical_Products.fasta "
    f"-o {str(AnnotaBasename)}_interproscan_hypothetical_output.gff3 "
    f"-f GFF3 -t p -goterms -iprlookup"
)

# Optionals
for variable in interpro:
    if str(interpro.get(variable)).lower() == "flag":
        interpro_command_line += f" -{str(variable)}"
    else:
        interpro_command_line += f" -{str(variable)} {str(interpro.get(variable))}"

logger.debug(str(interpro_command_line))

subprocess.getoutput(interpro_command_line)

# INTERPROSCAN parser (info_parser.py) can run without this result, but must be a valid file.
if os.path.isfile(str(AnnotaBasename + "_interproscan_hypothetical_output.gff3")) == 0:
    # Generate valid file
    open(f"{str(AnnotaBasename)}_interproscan_hypothetical_output.gff3", "w").close()
    logger.warning("INTERPROSCAN analysis return no results, moving on without this results.")
    logger.warning("Check if your sequences have special characters (like *), remove it and rerun")

logger.info("INTERPROSCAN finished for Hypothetical Proteins")
logger.info("Preparing file for INTERPROSCAN Annotated Proteins execution")

annotated_file = str(blast_folder / str(AnnotaBasename + "_annotated_products.txt"))

# Using only IDs from the file
annotated_id = [line.strip().split()[0] for line in open(annotated_file, "r")]
fasta_fetcher(str(augustus_folder / str("Clear_" + aug_parsing)), annotated_id,
              "Annotated_Products.fasta")

# Check if expected file exists
check_file("Annotated_Products.fasta")

logger.info("INTERPROSCAN Annotated Proteins file preparation complete")

logger.info("Running INTERPROSCAN with Annotated Proteins")

# INTERPROSCAN: commandline

# General
interpro_command_line = (
    f"{str(AnnotaPipeline.get('interpro_exe'))} -i Annotated_Products.fasta "
    f"-o {str(AnnotaBasename)}_interproscan_annotated_output.gff3 "
    f"-f GFF3 -t p -goterms -iprlookup"
)

for variable in config['INTERPROSCAN']:
    if str(interpro.get(variable)).lower() == "flag":
        interpro_command_line += f" -{str(variable)}"
    else:
        interpro_command_line += f" -{str(variable)} {str(interpro.get(variable))}"

logger.debug(str(interpro_command_line))

subprocess.getoutput(interpro_command_line)

# INTERPROSCAN parser (info_parser.py) can run without this result, but must be a valid file.
if os.path.isfile(str(AnnotaBasename + "_interproscan_annotated_output.gff3")) == 0:
    # Generate valid file
    open(f"{str(AnnotaBasename)}_interproscan_annotated_output.gff3", "w").close()
    logger.warning("INTERPROSCAN analysis return no results, moving on without this results.")
    logger.warning("Check if your sequences have special characters (like *), remove it and rerun")

logger.info("INTERPROSCAN finished for Annotated Proteins")

# HMMER ------------------------------------------------------------------------
logger = logging.getLogger('HMMSCAN')
logger.info("Running HMMSCAN with Hypothetical Proteins")

# General
hmmscan_command_line = (
    f"{str(AnnotaPipeline.get('hmm_exe'))} "
    f"--cpu {str(AnnotaPipeline.get('threads'))} "
    f"--tblout {str(AnnotaBasename)}_hmmscan_output.txt "
    f"--noali"
)

# Optionals
for variable in hmmscan:
    if str(hmmscan.get(variable)).lower() == "flag":
        hmmscan_command_line += f" --{str(variable)}"
    else:
        # These specific arguments are passed through '-'
        if any(arg == str(variable) for arg in ("E", "Z", "T")):
            hmmscan_command_line += f" -{str(variable)} {str(hmmscan.get(variable))}"
        # Everyone else are passed through '--'
        else:
            hmmscan_command_line += f" --{str(variable)} {str(hmmscan.get(variable))}"

hmmscan_command_line += (
    f" {str(AnnotaPipeline.get('pfam'))} "
    f"Hypothetical_Products.fasta "
    f"> /dev/null 2> hmmscan.err"
)

logger.debug(str(hmmscan_command_line))

subprocess.getoutput(hmmscan_command_line)

# Check if expected file exists
check_file(f"{str(AnnotaBasename)}_hmmscan_output.txt")

logger.info("HMMSCAN is finished")

# RPSBLAST ---------------------------------------------------------------------
logger = logging.getLogger('RPSBLAST')
logger.info("Running RPSBLAST with Hypothetical Proteins")

# General
rpsblast_command_line = (
    f"rpsblast -query Hypothetical_Products.fasta "
    f"-out {str(AnnotaBasename)}_rpsblast_output.outfmt6 "
    f"-db {str(AnnotaPipeline.get('cdd'))} "
    f"-outfmt \"6 qseqid sseqid sacc bitscore evalue ppos pident qcovs stitle\" "
    f"-num_threads {str(AnnotaPipeline.get('threads'))}"
)

# Optionals
for variable in rpsblast:
    if str(rpsblast.get(variable)).lower() == "flag":
        rpsblast_command_line += f" -{str(variable)}"
    else:
        rpsblast_command_line += f" -{str(variable)} {str(rpsblast.get(variable))}"

logger.debug(str(rpsblast_command_line))

subprocess.getoutput(rpsblast_command_line)

# Check if expected file exists
check_file(f"{str(AnnotaBasename)}_rpsblast_output.outfmt6")

logger.info("RPSBLAST is finished")

# -----------------------------------------------------------------------
logger = logging.getLogger('AnnotaPipeline')
logger.info("Parsing information from INTERPROSCAN, HMMSCAN and RPSBLAST")

subprocess.run([
    str(python_exe),
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
except Exception as warn:
    logger.warning("Failed to remove HMMSCAN log and temp dir")
    logger.debug(f"code error: {warn}")
    pass

logger.info("INTERPROSCAN, HMMSCAN and RPSBLAST execution and parsing is finished")

os.chdir(annota_pwd)

# -----------------------------------------------------------------------

logger.info("Preparing file for protein annotation")

subprocess.run([
    str(python_exe),
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
except Exception as warn:
    logger.warning("Failed to sort All_annotation_products.txt")
    logger.debug(f"code error: {warn}")
    pass

logger.info("All_Annotated_Products.txt file is complete")

# ------------  Defining what file will be used ---------------------------
if args.gff is not None and args.protein is not None:  # User gave protein file and gff file
    # Run parser to generate fasta_file
    gff_file = str(gff_path)
    logger.info(f"Generating AnnotaPipeline_{AnnotaBasename}_proteins.fasta")
    gfftofasta(str(python_exe))
    logger.info("Generating annotated GFF file")
    run_fastatogff(str(python_exe))
    logger.info(f"GFF file is ready - Check {AnnotaBasename}_Annotated_GFF.gff")
elif args.protein is not None and args.gff is None:  # User gave only protein file
    logger.info("GFF file wasn't given, FASTA file will only have annotations")
    logger.info("Generating simple FASTA")
    subprocess.run([
        str(python_exe),
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
    logger.info("GFF file wasn't given, skipping script fastatogff.py")
else:  # User selected run Augustus
    gff_file = augustus_folder / str("AUGUSTUS_" + str(AnnotaBasename) + ".gff")
    logger.info(f"Generating AnnotaPipeline_{AnnotaBasename}_proteins.fasta")
    gfftofasta(str(python_exe))
    logger.info("Generating annotated GFF file")
    run_fastatogff(str(python_exe))
    logger.info(f"GFF file is ready - Check {str(AnnotaBasename)}_Annotated_GFF.gff")

# Transfer annotation from proteins to Transcript (codingseq) file
annotate_codingseq(annota_pwd / str("AnnotaPipeline_" + AnnotaBasename + "_proteins.fasta"), 
    augustus_folder / str("AUGUSTUS_" + AnnotaBasename + ".codingseq"), AnnotaBasename)

logger.info("AnnotaPipeline has annotated the annotations on the annotated file.")

# -----------------------------------------------------------------------
# ---------------------- Kallisto ---------------------------------------
# Run kallisto
# If kalisto_method is empty, lack arguments for kallisto >> skip
# If proteins were given, lack cdscexon files >> skip
if kallisto_method == None or args.protein is not None:
    pass
else:
    kallisto_output_path = pathlib.Path(annota_pwd / str("4_TranscriptQuantification_" + AnnotaBasename))
    # kallisto_output_path = pathlib.Path(f"{annota_pwd / f'4_TranscriptQuantification_{AnnotaBasename}'}")
    # testei e funciona
    # documentacao do Pathlib
    # SubDeskTop = Path.joinpath(Desktop, "subdir")
    pathlib.Path(kallisto_output_path).mkdir(exist_ok=True)
    # Go to /4_TranscriptQuantification_
    os.chdir(kallisto_output_path)

    # Change logger
    logger = logging.getLogger('KALLISTO')

    # Annotated_Products.cdsexon 
    kallisto_run(str(python_exe),
        kallisto.get("kallisto_path"), kallisto_paired_end, kallisto_method,
        AnnotaBasename, f"{annota_pwd / f'AnnotaPipeline_{AnnotaBasename}_transcripts.fasta'}"
        f"-D{annota_pwd / f'AnnotaPipeline_{AnnotaBasename}_proteins.fasta'} " \
    )

    logger.info("KALLISTO execution and parsing is finished")

    # Return to AnnotaPipeline basedir
    os.chdir(annota_pwd)

# -----------------------------------------------------------------------
# ----------------------- Commet ----------------------------------------
if len(comet.get('comet_bash')) == 0:
    pass
else:
    # ----------------------------------- Path name  --------------------------------------------------
    if kallisto_method == None or args.protein is not None:
        comet_output_path = pathlib.Path(annota_pwd / str("4_PeptideIdentification" + AnnotaBasename))
    else:
        comet_output_path = pathlib.Path(annota_pwd / str("5_PeptideIdentification" + AnnotaBasename))
    # -------------------------------------------------------------------------------------------------
    # Change logger
    logger = logging.getLogger('COMET')
    # Go to /X_PeptideIdentification
    pathlib.Path(comet_output_path).mkdir(exist_ok=True)
    os.chdir(comet_output_path)
    # -------------------------------------------------------------------------------------------------
    # Check if overwrite parameters will be used
    if use_last_and_first == True:
        first_last_param = f"-F{comet.get('first')} -L{comet.get('last')}"
    else:
        first_last_param = str()
    # -------------------------------------------------------------------------------------------------
    mass_path = f"{str(comet.get('mass_files')).rstrip('/')}/"
    commet_command = f"{comet.get('comet_bash')} -P{comet.get('params')} " \
                     f"-D{annota_pwd / f'AnnotaPipeline_{AnnotaBasename}_proteins.fasta'} " \
                     f"{first_last_param} {mass_path}*"
    # -------------------------------------------------------------------------------------------------
    logger.info("COMET execution has started")
    logger.debug(commet_command)
    subprocess.getoutput(commet_command)
    logger.info("COMET execution is finished")
    logger.info("Parsing COMET output")
    # -------------------------------------------------------------------------------------------------
    # Get all output files from mass_path >> default output path
    file_names = pathlib.Path(mass_path).glob('*.pin')
    # Check_files for comet_output
    if len(file_names) == 0:
        logger.error("COMET returns no output")
        log_quit()
    # -------------------------------------------------------------------------------------------------
    logger = logging.getLogger('PERCOLATOR')
    logger.info("PERCOLATOR execution has started")
    # ----------------------------------PERCOLATOR RUN-------------------------------------------------
    for comet_output_file in file_names:
        # get only filename (without path), and remove comet range from filename (ex: filename.2-200.pin)
        percolator_out_basename = re.sub(r"\.[0-9].*","",comet_output_file.stem)
        percolator_command = f"{percolator.get('percolator_bash')} -r {percolator_out_basename}_peptide_output.tsv " \
                            f" -m {percolator_out_basename}_percolator_output.tsv" \
                            f" -B {percolator_out_basename}_decoy_output.tsv {comet_output_file}"
        logger.info(f"Runing percolator with sample: {comet_output_file}")
        logger.debug(percolator_command)
        subprocess.getoutput(percolator_command)

        logger.info(f"Parsing {percolator_out_basename}_percolator_output.tsv")
        check_file(f"{percolator_out_basename}_percolator_output.tsv")
        # --------------------------- PERCOLATOR PARSING -----------------------------------------------
        parser_percolator_command = f"{python_exe} {str(pipeline_pwd / 'percolator_parser.py')}" \
                            f" -p {percolator_out_basename}_percolator_output.tsv" \
                            f" -qv {percolator.get('qvalue')}" \
                            f" -b {AnnotaBasename}_{percolator_out_basename}"
        try:
            subprocess.getoutput(parser_percolator_command)
        except Exception as warn:
            logger.warning(f"Failed trying to parser {percolator_out_basename}_percolator_output.tsv")
            logger.debug(f"code error {warn}")
    logger.info("PERCOLATOR parsing is finished")
    # -------------------------------------------------------------------------------------------------
    logger.info("Creating quantitative report of Spectrum and Peptides")
    quantitative_proteomics(comet_output_path, AnnotaBasename)
    try:
        # Sort spectrum count
        os.system(f"sort -V {AnnotaBasename}_pre_total_Proteomics_Quantification.tsv " \
                  f"-o  {AnnotaBasename}_Total_Proteomics_Quantification.tsv")
        os.remove(f"{AnnotaBasename}_pre_total_Proteomics_Quantification.tsv")
    except Exception as warn:
        logger.warning("Failed to sort All_annotation_products.txt")
        logger.debug(f"code error: {warn}")
        pass

# Return to AnnotaPipeline basedir
os.chdir(annota_pwd)
# close logger
logging.shutdown()
