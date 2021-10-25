#!/usr/bin/python3

############################
###   DESCRIPTION HERE   ###
############################

# USAGE: python3 BLASTp_SwissProt_TriTrypDB.py \
#                    --seq Cleared_Linfantum_AUGUSTUS.aa \
#                    --swissprot uniprot_sprot.fasta \
#                    --tritrypdb TriTrypDB_proteins.fasta

"""---A AVENTURA VAI COMEÇAR-------------------------------------------------"""

import argparse
import subprocess
import logging
import sys
import re
from shutil import which

'''---ARGUMENTS AND [--help / -help / -h]------------------------------------'''

parser = argparse.ArgumentParser(
    add_help=False,  # removes original [--help]
    description='''Scritp to run and parse output from Swissprot and other database [NR or Trytrip],
    Please give at least one, through flags -spdb, -trbl or -nr    
    ''',
    epilog=""">>>> -nr, -spdb and -trbl are mutually exclusive arguments <<<<
    
Poof, you're a sandwich!""", formatter_class=argparse.RawTextHelpFormatter
)

requiredNamed = parser.add_argument_group('required arguments')
optionalNamed = parser.add_argument_group('optional arguments')

# mandatory arguments
#   type (default): string
requiredNamed.add_argument(
    '-s', '--seq', dest='seq',
    metavar='[protein_file]',
    help=('non-redundant protein fasta file'
          + ' to be used by the BLASTp suite'),
    required=True
)

requiredNamed.add_argument(
    '-sp', '--swissprot', dest='spdb',
    metavar='[UniProt_SwissProt_database]',
    help='destination to /SwissProt/database',
    required=True
)

requiredNamed.add_argument(
    '-basename', dest='basename',
    metavar='[It\'s a boy, and will be called Jonas]',
    help='basename',
    required=True
)
#   type (default): string
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument(
    '-trbl', '--Trembl', dest='trembl',
    metavar='[Trembl_database]',
    help='destination to Trembl database'
)

group.add_argument(
    '-spdb', '--specificdb', dest='specificdb',
    metavar='[SpecificDB_database]',
    help='destination to specific database, in EupathDB format',
)

group.add_argument(
    '-nr', '--nrdb', dest='nr',
    metavar='[NR_database]',
    help='destination to NR database',
)

optionalNamed.add_argument(
    '-id', '--identity', dest='id',
    metavar='', default=40,
    help=('Minimal identity to transfer annotation. If a Query have just results below this threshold '
          'will be considered as hypothetical'
          + ' (default: 40)')
)

optionalNamed.add_argument(
    '-cov', '--coverage', dest='cov',
    metavar='', default=30,
    help=('Minimal coverage to analise query result. Matchs below this threshold '
          'will not be considered'
          + ' (default: 30)')
)

optionalNamed.add_argument(
    '-kw', '--keywords', dest='keywords',
    metavar='[\"hypothetical,unspecified,fragment'
            ',partial,unknown,fragemnt\"]', default=str("hypothetical,unspecified,fragment,"
                                                        "partial,unknown,fragemnt"),
    help=('Keywords to search for hypothetical annotations. Please, pass each word followed by comma,'
            +' whitout spaces')
)

optionalNamed.add_argument(
    '-pos', '--positivity', dest='pos',
    metavar='', default=60,
    help=('Minimal positivity to transfer annotation. If a Query have just results below this threshold '
          'will be considered as hypothetical'
          + ' (default: 60)')
)

optionalNamed.add_argument(
    '-t', '--threads', dest='threads',
    metavar='', type=int, default=20,
    help='number of threads [int] (default: 20)'
)

optionalNamed.add_argument(
    '-hsps', dest='hsps',
    metavar='', type=int, default=10,
    help='max_hsps flag from blastp [int] (default: 10)'
)

optionalNamed.add_argument(
    '-evalue', dest='evalue',
    metavar='', type=float, default=0.00001,
    help='evalue flag from blastp [int] (default: 0.00001)'
)

# custom [--help] argument
optionalNamed.add_argument(
    '-h', '-help', '--help',
    action='help',
    default=argparse.SUPPRESS,  # hidden argument
    help='It\'s going to be legen - wait for it - dary!'
)

# arguments saved here
args = parser.parse_args()

# ----------------------- Create LogFile ------------------------------------

logger = logging.getLogger('Blast')

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
    filename="Blast.log",
    filemode='a'
)

stderr_logger = logging.getLogger('Blast')
sl = StreamToLogger(stderr_logger, logging.ERROR)
sys.stderr = sl

# --------------------------------------------------------------------------

'''---SOFTWARE SETTINGS------------------------------------------------------'''


class hit:

    def __init__(self, desc, bitscore):
        self.desc = desc
        self.bitscore = bitscore

    def __lt__(self, other):
        return self.bitscore > other.bitscore


def blast(arq1, arq2, db, hsps, evalue):
    datab = str(db)
    fmt = str("\"6 qseqid sseqid sacc bitscore"
              + " evalue ppos pident qcovs stitle\"")

    command = f"blastp -query {arq1} -out {arq2}" \
                    f" -db {datab} -evalue {evalue}" \
                    f" -outfmt {fmt}" \
                    f" -max_hsps {hsps}" \
                    f" -num_threads {str(args.threads)}"
    logger.info(command)
    subprocess.getoutput(command)


def temporary_query(arq):
    # temporary adding a line to the BLAST analysis
    #   so it accounts for the last query found by the BLAST algorithm
    temp = arq[-1].split("\t")
    del temp[0]
    temp.insert(0, "QueryTemp")
    arq.append("\t".join(temp))


def is_tool(name):
    if which(name) is None:
        logging.error(f"Program: {str(name)} must be avaliable in $PATH")
        logging.shutdown()
        sys.exit(1)


'''---Keywords-------------------------------------------------------------'''

# defining the keywords that will be used
#   to separate each HSP found in the BLAST output_file.txt:
keyword_list = args.keywords.split(",")


def parser_trembl(basename, result_blast, identidade, positividade, cov):
    ## Parser start

    swiss = open(str(result_blast), "r").read().splitlines()

    hyp = open(f"{str(basename)}_hypothetical_products.txt", "w")
    nhyp = open(f"{str(basename)}_annotated_products.txt", "a")
    all_anot = open(f"{str(basename)}_SpecifiedDB_annotations.txt", "w")
    # temporary adding a line
    #   so it accounts for the last query found by the BLAST analysis
    temporary_query(swiss)

    # creating the lists and counter that will be used
    #   to make sure the script goes through each HSP in every query
    # it's important that these lists are set back to NULL
    #   and the counter is set to zero before we start
    old_id = swiss[0].split("\t") 
    old_id = old_id[0]
    annots = []
    nhyp_list = []
    classification = []
    desc_list = []

    for query in swiss:
        line_split = query.split("\t")
        title = line_split[-1]  # get description
        title_split = title.split(" ", 1)[-1]  # get last part of description
        desc = title_split.split("OS=")[0].strip().rstrip()  # get only description
        new_id = line_split[0]
        # defining which file will receive each HSP depending on the counter number
        #   the script will write each HSP on its corresponding .txt file
        # IMPORTANT: the first time this loop runs it will add an empty line ("\n")
        #            to the first line of the hyp_file.txt
        if old_id != new_id:
            if "non_hypothetical" in ' '.join(classification):
                nhyp.write(f"{str(old_id)}\t")
                # Sort annotations by identity
                annots.sort()
                # Get best identity, first position of array
                nhyp.write(f"{str(annots[0].desc)}\n")
                nhyp.write(str('\n'))
                nhyp_list.append(str(old_id))

                all_anot.write(f"{old_id}\t{len(desc)} Annotation(s): [{';'.join(desc_list)}]\n")
            else:
                hyp.write(f"{str(old_id)}\n")
            # this just resets the count back to zero, before it starts again
            classification.clear()
            annots.clear()
            desc_list.clear()
        # if the HSP is a pseudogene (based on its title)
        #   it's counted separately from the other two
        else:
            if float(line_split[7]) > float(cov):
                # Check annotations, if any word doesn't match with keywords
                if not any(s in desc.lower() for s in keyword_list):
                    # Results that doesn't have match with keyword_list pass
                    # Check positivity and identity with subject
                    # If it's lower than threshold, it's not trustworthy, so, it's considered hypothetical
                    if float(line_split[5]) >= float(positividade) \
                            and float(line_split[6]) >= float(identidade):
                        # If annotation is strong, is considered as non_hypothetical
                        annots.append(hit(str(desc), float(line_split[3])))
                        desc_list.append(str(desc))
                        classification.append("non_hypothetical")
                    else:
                        classification.append("hypothetical")
                    # Check annotations, if any word match with keywords, is considered hypothetical
                else:
                    classification.append("hypothetical")
            else:
                pass
        # saving the information that will be written in the .txt files
        old_id = new_id

    hyp.close()
    nhyp.close()
    all_anot.close()


def parser_trytrip(basename, result_blast, identidade, positividade, cov):
    ## Parser start

    swiss = open(str(result_blast), "r").read().splitlines()

    hyp = open(f"{str(basename)}_hypothetical_products.txt", "w")
    nhyp = open(f"{str(basename)}_annotated_products.txt", "a")
    all_anot = open(f"{str(basename)}_SpecifiedDB_annotations.txt", "w")
    # temporary adding a line
    #   so it accounts for the last query found by the BLAST analysis
    temporary_query(swiss)

    # creating the lists and counter that will be used
    #   to make sure the script goes through each HSP in every query
    # it's important that these lists are set back to NULL
    #   and the counter is set to zero before we start
    old_id = swiss[0].split("\t")  # recebe a primeira query para começar a contagem
    old_id = old_id[0]
    annots = []
    nhyp_list = []
    classification = []
    desc_list = []

    for query in swiss:
        # split line by \t > separate columns > each line becomes a list
        title = query.split("\t")
        description = title[8].split("|")
        desc = description[5].replace("transcript_product=", "")
        new_id = title[0]
        # defining which file will receive each HSP depending on the counter number
        #   the script will write each HSP on its corresponding .txt file
        # IMPORTANT: the first time this loop runs it will add an empty line ("\n")
        #            to the first line of the hyp_file.txt
        if old_id != new_id:
            if "non_hypothetical" in ' '.join(classification):
                nhyp.write(f"{str(old_id)}\t")
                # Sort annotations by identity
                annots.sort()
                # Get best identity, first position of array
                nhyp.write(f"{str(annots[0].desc)}\n")
                nhyp_list.append(str(old_id))

                all_anot.write(f"{old_id}\t{len(desc)} Annotation(s): [{';'.join(desc_list)}]")
                all_anot.write(str('\n'))
            else:
                hyp.write(f"{str(old_id)}\n")
            # this just resets the count back to zero, before it starts again
            classification.clear()
            annots.clear()
            desc_list.clear()
        # if the HSP is a pseudogene (based on its title)
        #   it's counted separately from the other two
        else:
            if float(title[7]) > float(cov):
                # Check annotations, if any word doesn't match with keywords
                if not any(s in desc.lower() for s in keyword_list):
                    # Results that doesn't have match with keyword_list pass
                    # Check positivity and identity with subject
                    # If it's lower than threshold, it's not trustworthy, so, it's considered hypothetical
                    if float(title[5]) >= float(positividade) \
                            and float(title[6]) >= float(identidade):
                        # If annotation is strong, is considered as non_hypothetical
                        annots.append(hit(str(desc), float(title[3])))
                        desc_list.append(str(desc))
                        classification.append("non_hypothetical")
                    else:
                        classification.append("hypothetical")
                    # Check annotations, if any word match with keywords, is considered hypothetical
                else:
                    classification.append("hypothetical")
            else:
                pass
        # saving the information that will be written in the .txt files
        old_id = new_id

    hyp.close()
    nhyp.close()
    all_anot.close()


'''Parser nr, the same thing, but different'''


def parser_nr(basename, result_blast, identidade, positividade, cov):
    # Parser start

    nr = open(str(result_blast), "r").read().splitlines()

    hyp = open(f"{str(basename)}_hypothetical_products.txt", "w")
    nhyp = open(f"{str(basename)}_annotated_products.txt", "a")
    all_anot = open(f"{str(basename)}_SpecifiedDB_annotations.txt", "w")
    # temporary adding a line
    #   so it accounts for the last query found by the BLAST analysis
    temporary_query(nr)

    # creating the lists and counter that will be used
    #   to make sure the script goes through each HSP in every query
    # it's important that these lists are set back to NULL
    #   and the counter is set to zero before we start
    old_id = nr[0].split("\t")  # recebe a primeira query para começar a contagem
    old_id = old_id[0]
    annots = []
    classification = []
    desc_list = []

    for query in nr:
        # split line by \t > separate columns > each line becomes a list
        title = query.split("\t")
        description = title[8]
        # search for description without name in []
        desc = re.search(r'\s(.*?)\s\[.*', description).group(1)
        new_id = title[0]
        # defining which file will receive each HSP depending on the counter number
        #   the script will write each HSP on its corresponding .txt file
        # IMPORTANT: the first time this loop runs it will add an empty line ("\n")
        #            to the first line of the hyp_file.txt
        if old_id != new_id:
            if "non_hypothetical" in ' '.join(classification):
                nhyp.write(f"{str(old_id)}\t")
                # Sort annotations by identity
                annots.sort()
                # Get best identity, first position of array
                nhyp.write(f"{str(annots[0].desc)}\n")
                all_anot.write(f"{old_id}\t{len(desc)} Annotation(s): [{';'.join(desc_list)}]\n")
            else:
                hyp.write(f"{str(old_id)}\n")
            # this just resets the count back to zero, before it starts again
            classification.clear()
            annots.clear()
            desc_list.clear()
        # if the HSP is a pseudogene (based on its title)
        #   it's counted separately from the other two
        else:
            if float(title[7]) > float(cov):
                # Check annotations, if any word doesn't match with keywords
                if not any(s in desc.lower() for s in keyword_list):
                    # Results that doesn't have match with keyword_list pass
                    # Check positivity and identity with subject
                    # If it's lower than threshold, it's not trustworthy, so, it's considered hypothetical
                    if float(title[5]) >= float(positividade) \
                            and float(title[6]) >= float(identidade):
                        # If annotation is strong, is considered as non_hypothetical
                        classification.append("non_hypothetical")
                        annots.append(hit(str(desc), float(title[3])))
                        desc_list.append(desc)
                    else:
                        classification.append("hypothetical")
                # Check annotations, if any word match with keywords, is considered hypothetical
                else:
                    classification.append("hypothetical")
            else:
                pass
        # saving the information that will be written in the .txt files
        old_id = new_id

    hyp.close()
    nhyp.close()
    all_anot.close()


'''Parser swissprot, the same thing, but different'''


def process_swiss(basename, protein_seq, swiss_out, identidade, positividade, cov):
    # --------------------------Parser ----------------------------------------------------
    swiss = open(str(swiss_out), "r").read().splitlines()

    nhyp = open(f"{str(basename)}_annotated_products.txt", "w")
    swiss_anot = open(f"{str(basename)}_SwissProt_annotations.txt", "w")
    # temporary adding a line
    #   so it accounts for the last query found by the BLAST analysis
    temporary_query(swiss)

    # creating the lists and counter that will be used
    #   to make sure the script goes through each HSP in every query
    # it's important that these lists are set back to NULL
    #   and the counter is set to zero before we start
    old_id = swiss[0].split("\t")  # recebe a primeira query para começar a contagem
    old_id = old_id[0]
    nhyp_list = []
    classification = []
    desc = []
    annots = []

    # loop starting and parsing preparation
    for query in swiss:
        # split line by \t > separate columns > each line becomes a list
        line_split = query.split("\t")
        title = line_split[-1]  # get description
        title_split = title.split(" ", 1)[-1]  # get last part of description
        description = title_split.split("OS=")[0].strip().rstrip()  # get only description
        new_id = line_split[0]
        # defining which file will receive each HSP depending on the counter number
        #   the script will write each HSP on its corresponding .txt file
        # IMPORTANT: the first time this loop runs it will add an empty line ("\n")
        #            to the first line of the hyp_file.txt
        if old_id != new_id:
            if "non_hypothetical" in ' '.join(classification):
                nhyp.write(f"{str(old_id)}\t")
                nhyp_list.append(str(old_id))
                # Sort annotations by identity
                annots.sort()
                # Get best identity, first position of array
                nhyp.write(f"{str(annots[0].desc)}\n")

                nhyp_list.append(str(old_id))

                swiss_anot.write(f"{old_id}\t{len(desc)} Annotation(s): [{';'.join(desc)}]\n")

            desc.clear()
            classification.clear()
            annots.clear()
        else:
            if float(line_split[7]) > float(cov):
                # ifcov # mesmo q colocar flag qndo roda o blast -- apenas resultados acima rodar 30-90
                if not any(s in description.lower() for s in keyword_list):
                    if float(line_split[5]) >= float(positividade) and float(line_split[6]) >= float(identidade):
                        # If annotation is strong, is considered as non_hypothetical
                        classification.append("non_hypothetical")
                        annots.append(hit(str(description), float(line_split[3])))
                        desc.append(description)
                    else:
                        classification.append("hypothetical")
                else:
                    classification.append("hypothetical")
            else:
                pass
        # saving the information that will be written in the .txt files
        old_id = new_id

    # remove HSPs found by SwissProt from the original fasta_file input
    fasta = open(str(protein_seq), "r").read().split(">")

    for id_list in nhyp_list:
        for seq in fasta:
            if id_list in seq:
                fasta.remove(seq)

    new_fasta = open(f"{str(basename)}_BLASTp_AA_SwissProted.fasta", "w")
    new_fasta.write(">".join(fasta))
    new_fasta.close()
    swiss_anot.close()
    nhyp.close()


def no_hit(basename, blast6):

    # =============================== Parser sequences with no hit =============================
    # Get hit headers
    list_hit = set([line.strip().split()[0] for line in open(blast6, "r")])

    # Get all headers
    list_all = [
        line.strip().replace(">", "")
        for line in open(f"{str(basename)}_BLASTp_AA_SwissProted.fasta", "r")
        if line.startswith(">")
    ]

    no_hit_file = open(f"{str(basename)}_no_hit_products.txt", "w")
    for annotated in list_hit:
        if annotated in list_all:
            list_all.remove(annotated)
    if len(list_all) > 0:
        no_hit_file.write(f"{'\n'.join(list_all)}\n")
        no_hit_file.close()

    # =========================================================================================


def swiss_run():
    logger.info("Running BLAST against SwissProt")
    blast(args.seq, swiss_out, args.spdb, args.hsps, args.evalue)
    logger.info("Running parser SwissProt")


# Check BLAST, run Swissprot and parser it's results
is_tool("blastp")

# Run BLAST against swissprotDB
swiss_out = f"{str(args.basename)}_BLASTp_AAvsSwissProt.outfmt6"
swiss_run()
process_swiss(args.basename, args.seq, swiss_out, args.id, args.pos, args.cov)

# Secondary database
if args.nr is not None:
    odb_out_name = f"{str(args.basename)}_BLASTp_AAvsNRDB.outfmt6"
    logger.info("Running BLAST against NR")
    # Use the file above without sequences already annotated by swissprot
    blast(f"{args.basename}_BLASTp_AA_SwissProted.fasta", odb_out_name, args.nr, args.hsps, args.evalue)
    # ------------------------------
    logger.info("Running parser NR")
    parser_nr(args.basename, odb_out_name, args.id, args.pos, args.cov)
    logger.info("Parser NR done")
    # -------------No hit-----------
    no_hit(str(args.basename), odb_out_name)

elif args.trembl is not None:
    odb_out_name = f"{str(args.basename)}_BLASTp_AAvsTrembl.outfmt6"
    odb = args.trembl
    logger.info("Running BLAST against TrEMBL")
    # Use the file above without sequences already annotated by swissprot
    blast(f"{args.basename}_BLASTp_AA_SwissProted.fasta", odb_out_name, args.trembl, args.hsps, args.evalue)
    # ------------------------------
    logger.info("Running parser TrEMBL")
    parser_trembl(args.basename, odb_out_name, args.id, args.pos, args.cov)
    logger.info("Parser TrEMBL done")
    # ----------No hit--------------
    no_hit(str(args.basename), odb_out_name)

# EupathDB
elif args.specificdb is not None:
    odb_out_name = f"{str(args.basename)}_BLASTp_AAvsSpecifiedDB.outfmt6"
    odb = args.specificdb
    logger.info("Running BLAST against specificDB")
    # Use the file above without sequences already annotated by swissprot
    blast(f"{args.basename}_BLASTp_AA_SwissProted.fasta", odb_out_name, args.specificdb, args.hsps, args.evalue)
    # ------------------------------
    logger.info("Running parser specificDB")
    parser_trytrip(args.basename, odb_out_name, args.id, args.pos, args.cov)
    logger.info("Parser specificDB done")
    # ----------No hit--------------
    no_hit(str(args.basename), odb_out_name)
else:
    logger.error("Can't find any secondary database")
# ------------------------------
