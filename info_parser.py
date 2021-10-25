#!/usr/bin/env python3
import argparse
import os
import os.path
from shutil import ExecError

parser = argparse.ArgumentParser(
        add_help=False,  # removes original [--help]
        description='''Script to create annotations of hypothetical proteins with info from hmmer, interproscan and rpsblast

Input  - Sorted Terms [file_preparsed]
Ouput  - File with Id and annotation in the following way:

id      hypothetical protein (IPR[...],GO[...])

''',
        epilog="""Rise to fame, your time has come!""", formatter_class=argparse.RawTextHelpFormatter
)

requiredNamed = parser.add_argument_group('required arguments')
optionalNamed = parser.add_argument_group('optional arguments')

# mandatory arguments
#   type (default): string
requiredNamed.add_argument(
        '-ipr1', '--ipr1', dest='ipr1',
        metavar='[InterProScan_Out.txt]',
        help='Output Interpro',
        required=True
)

requiredNamed.add_argument(
        '-ipr2', '--ipr2', dest='ipr2',
        metavar='[InterProScan_Out.txt]',
        help='Output Interpro Info',
        required=True
)

requiredNamed.add_argument(
        '-a', '--annotated', dest='annot',
        metavar='[Annotated_products.txt]',
        help='File with id and annotations',
        required=True
)

requiredNamed.add_argument(
        '-nh', '--no_hit', dest='nohit',
        metavar='[No_hits.txt]',
        help='File with id of sequences with no hit in blast',
        required=True
)

requiredNamed.add_argument(
        '-hy', '--hypothetical', dest='hypo',
        metavar='[hypothetical_products.txt]',
        help='File with id of hypothetical proteins',
        required=True
)

# custom [--help] argument
optionalNamed.add_argument(
        '-h', '-help', '--help',
        action='help',
        default=argparse.SUPPRESS,  # hidden argument
        help='Ooh, Life is good... \
        As good as you wish!'
)

'''Output must be
id \t hypothetical protein (InterPro:IPR039726,GO:0004308,GO:0009405)
'''
# Input format fields
# id    DB      Acession        Annot   IPR     GO
# -------------------------------------------

# arguments saved here
args = parser.parse_args()


def temporary_query(arq):
        # temporary adding a line to the BLAST analysis
        #   so it accounts for the last query found by the BLAST algorithm
        temp = arq[-1].split("\t")
        del temp[0]
        temp.insert(0, "QueryTemp")
        arq.append("\t".join(temp))


# ------------------------------------------

def parser_interproscan(arq_entrada, arq_ipr):
        entrada = open(str(arq_entrada), "r").read().split("##FASTA")  # SPLITTING THE GFF3 FILE IN TWO CATEGORIES
        interp = entrada[0].split("##sequence-region")  # WE'LL BE USING ONLY THE FIRST PART OF THE GFF3 OUTPUT FILE
        del interp[0]
        ipr = open(str(arq_ipr), "a")

        lista = ["Coils", "MobiDBLite"] 

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
                        elif linha == seq_reg[1]:
                                # IGNORING THE FIRST LINE ON EACH GROUP OF QUERIES, AS IT'S NON-INFORMATIVE
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
                                try:
                                        evalue = str(linha[5]).replace(",", ".").lower()
                                except Exception:
                                        # Type of evalue format
                                        evalue = str(linha[5])
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
                                        ipr.write(f"{nome_subject}\t{db_certo}\t{name}\t{evalue}\t{anotacao_db}\t{interpro}\t{ontologia}\n")
        ipr.close()


# Function to write ids with no ipr and no GO (no result from InterproScan)
def write_no_ipr(ids_dict):
        # If there's any id without hit in interproscan file
        # It's will write here, without IPR or GO
        if len(ids_dict.keys()) > 0:
                for anot in ids_dict.keys():
                        output.write(f"{anot}\t{str(ids_dict.get(anot)).strip()}\n")
        if len(hypo) > 0:
                for hyp in hypo:
                        output.write(f"{hyp}\thypothetical protein\n")
        if len(nohit) > 0:
                for hyp in nohit:
                        output.write(f"{hyp}\thypothetical protein\n")
        output.close()


# Function to join IPRs and GOs from InterproScan with IDs hypothetical or without hits
def intepro_process(ids_dict):
        # File pre-precessed with IDs, IPRs and GOs
        interpro_out = open("Interpro_out_tmp.txt", "r").read().splitlines()
        temporary_query(interpro_out)

        old_id = interpro_out[0].split("\t")  # recebe a primeira query para começar a contagem
        old_id = old_id[0]

        # Initialize lists
        gos = []
        iprs = []
        for query in interpro_out:
                title = query.split("\t")
                new_id = title[0]
                ipr = title[5]
                go = title[6]
                if old_id != new_id:
                        # Sort IPRs and GOs
                        try:
                                iprs.sort(key=lambda item: item.split("IPR")[1])
                                gos.sort(key=lambda item: item.split("GO")[1])
                        except Exception:
                                # No GO or IP for this protein (not a true warning)
                                pass
                        # Check if interpro_result is in annotated
                        if old_id in ids_dict.keys():
                                if (len(iprs) > 0) and (len(gos) > 0):
                                        output.write(f"{old_id}\t{str(ids_dict.get(old_id)).strip()} ({str(','.join(iprs))},{str(','.join(gos))})\n")
                                elif (len(iprs) > 0) and (len(gos) == 0):
                                        output.write(f"{old_id}\t{str(ids_dict.get(old_id)).strip()} ({str(','.join(iprs))})\n")
                                elif (len(iprs) == 0) and (len(gos) > 0):
                                        output.write(f"{old_id}\t{str(ids_dict.get(old_id)).strip()} ({str(','.join(gos))})\n")
                                else:
                                        output.write(f"{old_id}\t{str(ids_dict.get(old_id)).strip()}\n")
                                del ids_dict[old_id]
                        # Else, interpro_result must be in hyphotetical
                        elif old_id in hypo:
                                if (len(iprs) > 0) and (len(gos) > 0):
                                        output.write(f"{old_id}\thypothetical protein ({str(','.join(iprs))},{str(','.join(gos))})\n")
                                elif (len(iprs) > 0) and (len(gos) == 0):
                                        output.write(f"{old_id}\thypothetical protein ({str(','.join(iprs))})\n")
                                elif (len(iprs) == 0) and (len(gos) > 0):
                                        output.write(f"{old_id}\thypothetical protein ({str(','.join(gos))})\n")
                                else:
                                        output.write(f"{old_id}\thypothetical protein\n")
                                hypo.remove(old_id)
                        else:
                                if (len(iprs) > 0) and (len(gos) > 0):
                                        output.write(f"{old_id}\thypothetical protein ({str(','.join(iprs))},{str(','.join(gos))})\n")
                                elif (len(iprs) > 0) and (len(gos) == 0):
                                        output.write(f"{old_id}\thypothetical protein ({str(','.join(iprs))})\n")
                                elif (len(iprs) == 0) and (len(gos) > 0):
                                        output.write(f"{old_id}\thypothetical protein ({str(','.join(gos))})\n")
                                else:
                                        output.write(f"{old_id}\thypothetical protein\n")
                                nohit.remove(old_id)
                        gos.clear()
                        iprs.clear()
                else:
                        if (ipr != "None") and (ipr not in iprs):
                                iprs.append(ipr)
                        if (go != "None") and (go not in gos):
                                gos.append(go)
                old_id = new_id

        write_no_ipr()
        os.remove("Interpro_out_tmp.txt")


# ---------------------- Create organized file with iprs and gos -----------------------
if (os.path.getsize(args.ipr1) == 0) and (os.path.getsize(args.ipr2) == 0):
        # Can't process empty file
        pass
else:
        # Almost one file has result - but, if one of them hasn't parser will crash
        if os.path.getsize(args.ipr1) != 0:
                parser_interproscan(args.ipr1, "Interpro_out_tmp.txt")
        if os.path.getsize(args.ipr2) != 0:
                parser_interproscan(args.ipr2, "Interpro_out_tmp.txt")

# ---------------------- Pre-parse annotated products -------------------------
annot = open(args.annot, "r").read().splitlines()
ids = []
desc = []
for anot in annot:
        anot = anot.split("\t")
        ids.append(anot[0])
        desc.append(anot[1])

ids_dict = dict(zip(ids, desc))  # Create dictionary with id as key and annotation as value

# ---------------------- Create list with hypothetical ids -----------------
if os.path.getsize(args.hypo) != 0:
        hypo = open(args.hypo, "r").read().splitlines()
else:
        hypo = []
# ---------------------- Create list with no_hits ids -----------------
if os.path.getsize(args.nohit) != 0:
        nohit = open(args.nohit, "r").read().splitlines()
else:
        nohit = []
# Create output
output = open("All_annotation_products.txt", "w")

if os.path.getsize(args.ipr1) == 0 and os.path.getsize(args.ipr2) == 0:
        write_no_ipr(ids_dict)
else:
        intepro_process(ids_dict)
