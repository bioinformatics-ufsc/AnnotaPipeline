#!/usr/bin/env python3
import argparse
import os
import os.path

# parser = argparse.ArgumentParser(
#         add_help=False,  # removes original [--help]
#         description='''Script join information from Annotations, Functional annotations, transcriptomics and Proteomics

# Expected output will be a xlsx file like with this columns:

# Protein_ID  Annotation  Superfamily IPR GO  Transcript (present|absent) Total_peptide Unique_Peptide

# ''',
#         epilog="""Rise to fame, your time has come!""", formatter_class=argparse.RawTextHelpFormatter
# )

# requiredNamed = parser.add_argument_group('required arguments')
# optionalNamed = parser.add_argument_group('optional arguments')

# # mandatory arguments
# #   type (default): string
# requiredNamed.add_argument(
#         '-annot', dest='annot',
#         metavar='[All_annotation.txt]',
#         help='Annotapipeline Output containing proteins annotated',
#         required=True
# )

# requiredNamed.add_argument(
#         '-ipr_annot', dest='ipr_annot',
#         metavar='[interproscan_annotated_output.gff3]',
#         help='Output from InterproScan with Annotated Proteins',
#         required=True
# )

# requiredNamed.add_argument(
#         '-ipr_hyp', dest='ipr_hyp',
#         metavar='[interproscan_hypothetical_output.gff3]',
#         help='Output from InterproScan with hypothetical Proteins',
#         required=True
# )

# optionalNamed.add_argument(
#         '-tr', dest='tr',
#         metavar='[interproscan_Transcript_Quantification.tsv]',
#         help='Summary from transcriptomics analysis (optional)',
# )

# optionalNamed.add_argument(
#         '-proteomics', dest='proteomics',
#         metavar='[interproscan_Total_Proteomics_Quantification.tsv]',
#         help='Summary from proteomics analysis (optional)',
# )

# # custom [--help] argument
# optionalNamed.add_argument(
#         '-h', '-help', '--help',
#         action='help',
#         default=argparse.SUPPRESS,  # hidden argument
#         help='Ooh, Life is good... \
#         As good as you wish!'
# )

# # arguments saved here
# args = parser.parse_args()


superfamily_dict = {}
ipr_dict = {}
go_dict = {}
annotation_dict = {}


def get_interpro_info(arq_entrada):
        entrada = open(str(arq_entrada), "r").read().split("##FASTA")  # SPLITTING THE GFF3 FILE IN TWO CATEGORIES
        interp = entrada[0].split("##sequence-region")  # WE'LL BE USING ONLY THE FIRST PART OF THE GFF3 OUTPUT FILE
        del interp[0]

        unwanted_db = ["Coils", "MobiDBLite"] 

        def create_or_reset_lists():
                ontologia = []
                interpro = []
                superfamily = []
                return ontologia, interpro, superfamily

        def add_to_dicts(id):
                superfamily_dict[id] = ",".join(superfamily)
                ipr_dict[id] = ",".join(interpro)
                go_dict[id] = ",".join(ontologia)

        def add_none_to_empty():
                if not ontologia:
                        ontologia.append("None")
                if not interpro:
                        interpro.append("None")
                if not superfamily:
                        interpro.append("None")

        # TREATING EACH QUERY, RETRIEVING THE INFORMATION THAT WILL BE WRITTEN ON THE FIRST OUTPUT FILE
        for seq_reg in interp:
                ontologia, interpro, superfamily = create_or_reset_lists()
                seq_reg = seq_reg.splitlines()
                for linha in seq_reg:
                        if linha == seq_reg[0]:
                                linha = linha.split(" ")
                                del linha[0]
                        elif linha == seq_reg[1]:
                                if linha == seq_reg[-1]:
                                        add_none_to_empty()
                                        add_to_dicts(nome_subject)
                                        ontologia, interpro, superfamily = create_or_reset_lists()
                                # IGNORING THE FIRST LINE ON EACH GROUP OF QUERIES, AS IT'S NON-INFORMATIVE
                                pass
                        else:
                                linha = linha.split("\t")
                                nome_subject = linha[0]
                                db = linha[1]
                                if not any(s in db for s in unwanted_db):
                                        anotation = linha[-1].split(";")
                                        for field in anotation:
                                                if "Ontology" in field:
                                                        ontologia.append(field.replace('"', "").replace("Ontology_term=", ""))
                                                if "Dbxref" in field:
                                                        interpro.append(field.replace('"', "").replace("Dbxref=", ""))
                                        if "SUPERFAMILY" in db:
                                                superfamily_line = linha[-1].split(";")
                                                for field in superfamily_line:
                                                        if "Name" in field:
                                                                superfamily.append(field.replace('"', "").replace("Name=", ""))
                                if '\t'.join(linha) == seq_reg[-1]:
                                        add_none_to_empty()
                                        add_to_dicts(nome_subject)
                                        ontologia, interpro, superfamily = create_or_reset_lists()


get_interpro_info("TrangeliSC58_interproscan_hypothetical_output.gff3")
#get_interpro_info("TrangeliSC58_interproscan_annotated_output.gff3")

