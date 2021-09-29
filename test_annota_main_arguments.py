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
import re
import pathlib
import shutil
import subprocess
import sys

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



def log_quit():
    logger.info("Exiting")
    logging.shutdown()
    sys.exit(1)


def annotate_codingseq(aa_fasta, codingseq_fasta):
    # Linha de anotacao completa
    anno_all = [line.strip() for line in open(aa_fasta) if ">" in line]
    anno_all = [sub.replace(">", "") for sub in anno_all]
    # IDs simplificadas
    anno_ids = [name.split("|")[0].strip() for name in anno_all]

    corrrected_ids = dict(zip(anno_ids, anno_all))

    id_dict  = SeqIO.to_dict(SeqIO.parse(codingseq_fasta, "fasta"))

    corrected_fasta = str("celegans_annotated_seqs.fasta")

    with open(corrected_fasta, "w") as corrected:
        for key, record in id_dict.items():
            for id_key in corrrected_ids.keys():
                if id_key in key:
                    record.id = corrrected_ids.get(id_key)
                    record.description = ""
                    SeqIO.write(record, corrected, "fasta")


def kallisto_run(kallisto_exe, paired_end, method, basename, fasta, output_type):
    
    kallisto_command_index = f"{kallisto_exe} index -i {basename}_kallisto_index_{output_type}.idx {fasta}"
    logger.info(f"Running Kallisto index with {output_type}")
    logger.info(f"{kallisto_command_index}")
    subprocess.getoutput(kallisto_command_index)

    # Standart command line for kallisto index
    # kallisto index -i transcripts.idx transcripts.fasta
    if paired_end == True:
        if len(kallisto.get("l")) != 0:
            l_flag = f"-l {kallisto.get('l')}"
            logger.info(f"Kallisto will run with -l {kallisto.get('l')}")
        else:
            l_flag = ""
        if len(kallisto.get("s")) != 0:
            s_flag = f"-s {kallisto.get('s')}"
            logger.info(f"Kallisto will run with -s {kallisto.get('s')}")
        else:
            s_flag = ""
       
        kallisto_command_quant = (
            f"{kallisto_exe} quant -i {basename}_kallisto_index_{output_type}.idx "
            f"{s_flag} {l_flag} "
            f"-o kallisto_output_{output_type} -b {kallisto.get('bootstrap')} {kallisto.get('rnaseq-data')}"
        )
        logger.info(f"Running Kallisto quant with {output_type}")
        logger.info(kallisto_command_quant)
        subprocess.getoutput(kallisto_command_quant)
        # Standart command line for kallisto quant paired end
        # kallisto quant -i transcripts.idx -o output -b 100 reads_1.fastq reads_2.fastq
    else:
        kallisto_command_quant = (
            f"{kallisto_exe} quant -i {basename}_kallisto_index_{output_type}.idx "
            f"-l {kallisto.get('s')} -s {kallisto.get('l')} --single "
            f"-o kallisto_output_{output_type} -b {kallisto.get('bootstrap')} {kallisto.get('rnaseq-data')}"
        )
        logger.info(f"Running Kallisto quant with {output_type}")
        logger.info(kallisto_command_quant)
        subprocess.getoutput(kallisto_command_quant)
        # Standart command line for kallisto quant single end
        # kallisto quant -i transcripts.idx -o output -b 100 --single -l 180 -s 20 reads_1.fastq

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
        f"python3 {kallisto_parser_path} "
        f"-ktfile kallisto_output/abundance.tsv "
        f"-basename {AnnotaBasename}_{output_type} {kallisto_parser_flag}"
    )
    logger.info(f"Running Parser for Kallisto with {output_type}")
    logger.info(kallisto_parser_command)
    subprocess.getoutput(kallisto_parser_command)


def kallisto_check_parameters():
    # kallisto method indicate wich method to parse (mean, median, value) 
    # kallisto paired_end indicate rnaseq type (important to run kallisto)
    global kallisto_method, kallisto_paired_end
    # This variable will store method to parse kallisto output
    kallisto_method = None
    # if kallisto path were given, check other arguments, if not pass.
    if len(config[str('KALLISTO')].get("kallisto_path")) == 0:
        logger.info("Arguments for Kallisto are empty. This step will be skipped.")
    else:
        kallisto_check = []
        for argument in ("rnaseq-data", "median", "mean", "value"):
            if argument == "rnaseq-data":
                if len(config[str('KALLISTO')].get("rnaseq-data").split()) > 2:
                    logger.error(
                        "Error, there are more arguments than required for rnaseq-data (KALLISTO). " \
                        "Pass one if your data is from single-end, " \
                        "and two files if your data is from paired end!"
                    )
                    log_quit()
                ###### This box check how many files were given in rnaseq-data #####
                elif len(config['KALLISTO'].get("rnaseq-data").split()) == 2:
                    kallisto_paired_end = True
                    logger.info(f"Kallisto will run with paired end data")
                    kallisto_check.append(argument)
                elif len(config['KALLISTO'].get("rnaseq-data").split()) == 1:
                    kallisto_paired_end = False
                    # Check required arguments for single end data
                    if len(config['KALLISTO'].get('l')) == 0:
                        logger.error("Error. Mandatory argument for single end data 'l' is empty")
                        log_quit()
                    elif len(config['KALLISTO'].get('s')) == 0:
                        logger.error("Error. Mandatory argument for single end data 's' is empty")
                        log_quit()
                    logger.info(f"Kallisto will run with single end data")
                    kallisto_check.append(argument)
                else:
                    logger.error(
                        "Error, check values for rnaseq-data (KALLISTO). " \
                        "Pass one if your data is from single-end, " \
                        "and two files if your data is from paired end!"
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
            logger.info("Arguments for Kallisto are empty. This step will be skipped.")
            kallisto_method = None
        # Rna-seq data is required for calisto
        if 'rnaseq-data' in kallisto_check:
            # if there is rna-seq data, check if method is correctly given
            if len(kallisto_check) > 2:
                logger.error("Error, there is more than one method selected to parse kallisto ouput. Please, review .config file.")
                log_quit()
            elif len(config[str('KALLISTO')].get("bootstrap")) == 0:
                logger.error("Kallisto bootstrap is empty, default value is 0. At least pass this value")
                log_quit()
            else:
                logger.info(f"Kallisto will run with method: {kallisto_check[1]}")
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
        if len(config['COMET'].get('last'))!=0:
            last_check = True
        # ----------------------------------------------

        for argument in ("params", "mass_files"):
            if len(config['COMET'].get(argument)) == 0:
                logger.error(f"Parameter [{argument}] from section [COMET] is null")
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
            logger.info("WARNING [COMET]: values for 'fist' and 'last' will overwrite those in params file")
            global use_last_and_first
            use_last_and_first = True
        else:
            use_last_and_first = False
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
    if sum([sp_verify, nr_verify, trembl_verify]) == 0:
        logger.error("Error: there is no secondary database. Please review config file!")
        log_quit()
    if sum([sp_verify, nr_verify, trembl_verify]) == 2:
        logger.error("Error: there are two secondary databases. Select one of them in the config file.")
        log_quit()


def fasta_fetcher(input_fasta, id_list, fetcher_output):
    wanted = sorted(set(id_list))
    records = (seq for seq in SeqIO.parse(input_fasta, "fasta") for r in wanted if r in seq.id)
    count = SeqIO.write(records, fetcher_output, "fasta")
    if count < len(wanted):
        logger.warning("IDs not found in input FASTA file")


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
        logger.info("Sometimes things don't work as we expect... and that's ok. "
                    "But, seriously, check your inputs and outputs and rerun.")


def is_tool(name):
    if which(name) is None:
        logging.error("Program: " + str(name) + " must be avaliable in $PATH")
        logging.shutdown()
        sys.exit(1)


def add_features(feature, data, save):
    parcial_feature = data.groupby(['ProteinID', f"{feature}"]).size().sort_values(ascending=False).reset_index(name=f'Unique {feature}')
    parcial_feature_2 = parcial_feature.groupby(['ProteinID']).size().sort_values(ascending=False).reset_index(name=f'Total {feature}')
    feature_process = parcial_feature.loc[parcial_feature[f'Unique {feature}'] == 1].drop(columns=[f"{feature}"])
    feature_process_2 = feature_process.groupby("ProteinID").count().reset_index()
    if 'ProteinID' not in save.columns:
        save["ProteinID"] = parcial_feature_2["ProteinID"]
    # Add unique feature column
    save = save.set_index("ProteinID").join(feature_process_2.set_index("ProteinID")).reset_index()
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
    total.to_csv(f"{basename}_Total_Proteomics_Quantification.tsv", sep="\t", index=False)


# --- CHECK EACH BOX OF VARIABLES ----------------------------------------------

sections_config = config.sections()
check_parameters(sections_config)

# --- SCRIPT LOCATION ----------------------------------------------------------

script_pwd = pathlib.Path(sys.argv[0]).absolute()

# AnnotaPipeline location
pipeline_pwd = script_pwd.parents[0]

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

# initial directory
home_dir_pwd = pathlib.Path.cwd()

home_dir = f"AnnotaPipeline_{str(AnnotaBasename)}"

pathlib.Path(home_dir).mkdir(exist_ok=True)

annota_pwd = pathlib.Path(home_dir_pwd / home_dir)

pathlib.Path(annota_pwd).mkdir(exist_ok=True)

#-----------------------------------------------------------
# ----------------------- Commet ----------------------------------------
if len(comet.get('comet_bash')) == 0:
    pass
else:
    if kallisto_method == None or args.protein is not None:
        comet_output_path = pathlib.Path(annota_pwd / str("4_PeptideIdentification_" + AnnotaBasename))
    else:
        comet_output_path = pathlib.Path(annota_pwd / str("5_PeptideIdentification_" + AnnotaBasename))
    logger = logging.getLogger('COMET')
    # Mudar loogger no log, facilita identificacao pra debug
    # Go to /X_PeptideIdentification
    pathlib.Path(comet_output_path).mkdir(exist_ok=True)
    os.chdir(comet_output_path)

    # Check if overwrite parameters will be used
    if use_last_and_first == True:
        first_last_param = f"-F{comet.get('first')} -L{comet.get('last')}"
    else:
        first_last_param = str()
    # Mass files location
    mass_path = f"{str(comet.get('mass_files')).rstrip('/')}/"
    commet_command = f"{comet.get('comet_bash')} -P{comet.get('params')} " \
                f"-D{annota_pwd / f'AnnotaPipeline_{AnnotaBasename}_proteins.fasta'} " \
                f"{first_last_param} {mass_path}/*.{comet.get('mass_files_ext')}"

    logger.info("COMET execution has started")
    logger.info(commet_command)
    subprocess.getoutput(commet_command)
    logger.info("COMET execution is finished")
    # ----------- Create Comet Output path --------------------------
    comet_path = pathlib.Path(comet_output_path / str("COMET_Output"))
    pathlib.Path(comet_path).mkdir(exist_ok=True)
    # ----------------------------------------------------------------
    # Get all output files from mass_path >> default output path
    file_names = pathlib.Path(mass_path).glob('*.pin')
    # ----------- Create Percolator Output path ----------------------
    logger = logging.getLogger('PERCOLATOR')
    logger.info("PERCOLATOR execution has started")
    # ----------------------------------------------------------------
    # ----------- Create Percolator Output path ----------------------
    percolator_path_raw = pathlib.Path(comet_output_path / str("Percolator_RAW"))
    pathlib.Path(percolator_path_raw).mkdir(exist_ok=True)
    # ----------------------------------------------------------------
    percolator_path_parsed = pathlib.Path(comet_output_path / str("Percolator_PARSED"))
    pathlib.Path(percolator_path_parsed).mkdir(exist_ok=True)
    # ----------------------------------------------------------------

    for comet_output_file in file_names:
        # -----------------------------------------
        # Copy Comet output files
        shutil.copy(comet_output_file, comet_path)
        # -----------------------------------------
        # --------- RUN Percolator inside Percolator RAW path -------------------
        os.chdir(percolator_path_raw)
        percolator_out_basename = re.sub(r"\.[0-9].*","",comet_output_file.stem)
        percolator_command = f"{percolator.get('percolator_bash')} -r {percolator_out_basename}_peptide_output.tsv" \
                             f" -m {percolator_out_basename}_percolator_output.tsv" \
                             f" -B {percolator_out_basename}_decoy_output.tsv {comet_output_file}"
        logger.info(f"Runing percolator with sample: {comet_output_file}")
        subprocess.getoutput(percolator_command)
        logger.info(f"Parsing {percolator_out_basename}_percolator_output.tsv")
        check_file(f"{percolator_out_basename}_percolator_output.tsv")
        # --------- RUN Percolator parser inside Percolator PARSED path -------------------
        os.chdir(percolator_path_parsed)
        parser_percolator_command = f"{python_exe} {str(pipeline_pwd / 'percolator_parser.py')}" \
                            f" -p {percolator_path_raw}/{percolator_out_basename}_percolator_output.tsv" \
                            f" -qv {percolator.get('qvalue')}" \
                            f" -b {AnnotaBasename}_{percolator_out_basename}"
        try:
            subprocess.getoutput(parser_percolator_command)
        except Exception as warn:
            logger.warning(f"Failed trying to parser {percolator_out_basename}_percolator_output.tsv")
            logger.debug(f"code error {warn}")
        os.chdir(comet_output_path)

    logger.info("PERCOLATOR parsing is finished")
    logger.info("Creating quantitative report of Spectrum and Peptides")
    quantitative_proteomics(f"{percolator_path_parsed}", AnnotaBasename)

# Return to AnnotaPipeline basedir
os.chdir(annota_pwd)
# close logger
logging.shutdown()