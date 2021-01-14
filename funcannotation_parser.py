# _*_ coding: utf-8 _*_

import re
import logging
import os
import argparse
import sys

# ---------------Parser arguments ----------------
parser = argparse.ArgumentParser(
    add_help=False,  # removes original [--help]
    description='''    
    Script to parser Interproscan, RPSblast and Pfam results.
    WARNING: Results from Coils and MobiDBLite won't be parsed
    ''',
    epilog="""And shall the hopeful words bring love inside your heart ...""",
    formatter_class=argparse.RawTextHelpFormatter
)

requiredNamed = parser.add_argument_group('required arguments')
optionalNamed = parser.add_argument_group('optional arguments')

# mandatory arguments
#   type (default): string
requiredNamed.add_argument(
    '-interpro', '--inter', dest='interpro',
    metavar='[InterProScan_Output.gff3]',
    help='InterProScan_Input',
    required=True
)

requiredNamed.add_argument(
    '-hmm', '--hmmscan', dest='hmm',
    metavar='[HMMSCAN_Output.txt]',
    help='HMMSCAN_Input',
    required=True
)

requiredNamed.add_argument(
    '-rpsblast', '--rpsblast', dest='rpsblast',
    metavar='[RPSblast_Output.txt]',
    help='RPSblast_Input',
    required=True
)

requiredNamed.add_argument(
    '-basename', '--basename', dest='basename',
    metavar='[It\'s a boy, and will be called Jonas]',
    help='basename',
    required=True
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

logger = logging.getLogger('Functional Annotation')

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
    filename="Functional_Annotation.log",
    filemode='a'
)

stderr_logger = logging.getLogger('Functional Annotation')
sl = StreamToLogger(stderr_logger, logging.ERROR)
sys.stderr = sl

# --------------------------------------------------------------------------


def parser_interproscan(arq_entrada, arq_ipr, arq_saida):
    entrada = open(str(arq_entrada), "r").read().split("##FASTA")  # SPLITTING THE GFF3 FILE IN TWO CATEGORIES
    interp = entrada[0].split("##sequence-region")  # WE'LL BE USING ONLY THE FIRST PART OF THE GFF3 OUTPUT FILE
    del interp[0]
    ipr = open(str(arq_ipr), "w")
    saida = open(str(arq_saida), "a")

    lista = ["Coils", "MobiDBLite"]  # THIS IS SPECIFIC TO MY ANALYSIS: I'M NOT LOOKING FOR STRUCTURAL EVIDENCE

    # TREATING EACH QUERY, RETRIEVING THE INFORMATION THAT WILL BE WRITTEN ON THE FIRST OUTPUT FILE
    for seq_reg in interp:
        seq_reg = seq_reg.splitlines()
        for linha in seq_reg:
            if linha == seq_reg[0]:
                linha = linha.split(" ")
                del linha[0]
                # nome_query = linha[0]
                # start_query = linha[1]
                # stop_query = linha[2]
            elif linha == seq_reg[1]:  # IGNORING THE FIRST LINE ON EACH GROUP OF QUERIES, AS IT'S NON-INFORMATIVE
                pass
            else:
                ontologia = str(None)
                anotacao_db = str(None)
                name = str(None)
                interpro = str(None)
                linha = linha.split("\t")
                # start_query = linha[3]
                # stop_query = linha[4]
                nome_subject = linha[0]
                db = linha[1]
                if not any(s in db for s in lista):
                    db_certo = db
                    anotacao = linha[-1].split(";")
                    for a in anotacao:
                        if "Ontology" in a:
                            ontologia = a.replace('"', "").replace("Ontology_term=", "")
                        if "signature_" in a:
                            anotacao_db = a.replace("signature_desc=", "")
                        if "Name" in a:
                            name = a.replace('"', "").replace("Name=", "")
                        if "Dbxref" in a:
                            interpro = a.replace('"', "").replace("Dbxref=", "")
                    ipr.write(str(nome_subject) + "\t" + str(db_certo) + "\t" + str(name) + "\t" +
                              str(anotacao_db) + "\t" + str(interpro) + "\t" + str(ontologia) + "\n")
                    saida.write(str(nome_subject) + "\t" + str(db_certo) + "\t" + str(name) + "\t" +
                                str(anotacao_db) + "\t" + str(interpro) + "\t" + str(ontologia) + "\n")


def pfam_format(arq_entrada, arq_saida):
    entrada = open(str(arq_entrada), "r").read().splitlines()
    saida = open(str(arq_saida), "w")
    del entrada[0:3]  # Del cabeçalho
    del entrada[-10:]  # Del rodapé
    for linha in entrada:
        linha = re.sub(r'\s+', ' ', linha).split(" ")  # retira todos os espaços e quebra as colunas
        aux = " ".join(linha[18:])  # reune as palavras por espaço
        del linha[18:]  # depois da coluna 18 apenas palavras
        linha.append(aux)  # reescreve a linha
        saida.write("\t".join(linha) + "\n")
    saida.close()


def parser_pfam(arq_entrada, arq_saida):
    entrada = open(str(arq_entrada), "r").read().splitlines()
    saida = open(str(arq_saida), "a")
    for linha in entrada:
        linha = linha.strip()
        linha = linha.split("\t")
        query = linha[2]
        db = "Pfam"  # ESPECIFICALLY WRITING THIS, AS THERE ARE NO OTHER DBs ON THIS PART OF THE ANALYSIS
        ontologia = str(None)
        anotacao_db = linha[-1]
        name = linha[1]
        interpro = str(None)
        saida.write(str(query) + "\t" + str(db) + "\t" + str(name) + "\t" + str(anotacao_db) + "\t" +
                    str(interpro) + "\t" + str(ontologia) + "\n")


def parser_rpsblast(arq_entrada, arq_rps, arq_saida):
    entrada = open(str(arq_entrada), "r").read().splitlines()
    rps = open(str(arq_rps), "w")
    saida = open(str(arq_saida), "a")
    for linha in entrada:
        hit = linha.split("\t")
        query = hit[0]
        db = "CDD"
        name = hit[1]  # Acesso_DB
        anotacao = hit[-1]
        interpro = str(None)
        ontologia = str(None)
        rps.write(
            str(query) + "\t" + str(db) + "\t" + str(name) + "\t" + str(anotacao) + "\t" + str(interpro) +
            "\t" + str(ontologia) + "\n")
        saida.write(
            str(query) + "\t" + str(db) + "\t" + str(name) + "\t" + str(anotacao) + "\t" + str(interpro) +
            "\t" + str(ontologia) + "\n")


def sort_arq(arq_entrada, arq_saida):
    hits = open(str(arq_entrada), "r")
    ordem = open(str(arq_saida), "w")
    ordem.write("ID\tDB\tDB_ACCESS\tDESCRIPTION\tIPR\tGO\n")
    lines = hits.readlines()
    lines.sort()
    for line in lines:
        ordem.write(str(line))


parser_interproscan(args.interpro, str("InterProScan_Out_" + args.basename + ".txt"),
                    str("Temp_" + args.basename + ".txt"))
logger.info("Interpro parser done")
pfam_format(args.hmm, str("Hmmscan_Out_" + args.basename + ".txt"))
logger.info("Hmmscan format done")
parser_pfam(str("Hmmscan_Out_" + args.basename + ".txt"), str("Temp_" + args.basename + ".txt"))
logger.info("Hmmscan parser done")
parser_rpsblast(args.rpsblast, str("RPSblast_Out_" + args.basename + ".txt"), str("Temp_" + args.basename + ".txt"))
logger.info("RPSblast parser done")
sort_arq(str("Temp_" + args.basename + ".txt"), str(args.basename + "_Grouped_Hypothetical_Information.txt"))
logger.info("Sorting querys in files")
os.system("rm " + str("Temp_" + args.basename + ".txt"))
logger.info("Temporary file removed")
logger.info("Script finished")
