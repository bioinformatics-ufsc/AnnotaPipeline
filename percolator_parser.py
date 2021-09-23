#!/usr/bin/python3

# --- A AVENTURA VAI COMEÇAR ---------------------------------------------------

import argparse
import pandas as pd
import sys
from pathlib import Path
import re

# PARSER ARGUMENTS -------------------------------------------------------------

parser = argparse.ArgumentParser(
    add_help=False,
    description='''
Script to parse Percolator Output
    ''',
    epilog="""Shadows will fade... someday...""",
    formatter_class=argparse.RawTextHelpFormatter
)

required_args = parser.add_argument_group("required arguments")
optional_args = parser.add_argument_group("custom arguments")

# mandatory arguments
required_args.add_argument(
    "-p", "--percolator-output", dest="proteom",
    type=argparse.FileType('r'), 
    metavar="[Percolator output]",
    help="Peptide quantification file (Percolator output)",
    required=True
)

required_args.add_argument(
    "-b", "--basename", dest="basename",
    metavar="[Jonas]",
    type=str, 
    help="It's a boy and will be called Jonas",
    required=True
)

# custom [--help] argument
required_args.add_argument(
    "-qv", "-q-value", dest='qvalue',
    metavar = "[between 0-1]",
    type=float, 
    help="Min qvalue to accept a peptide",
    required=False
)

# custom [--help] argument
optional_args.add_argument(
    "-h", "-help", "--help",
    action="help",
    default=argparse.SUPPRESS,
    help=""
)


# arguments saved here
args = parser.parse_args()

# Check qvalue
if  0 < args.qvalue > 1:
    print("-qv value must be between 0-1")
    sys.exit(1)

# --- PARSE QVALUE -----------------------------

df = pd.read_csv(args.proteom, sep="\t", skiprows=1, index_col=False)
# If needed, filter by charge


with pd.option_context('mode.chained_assignment', None):
    df_parsed = df.loc[df['q-value'] <= float(args.qvalue)]

# Save all
df_parsed.to_csv(f"{args.basename}_raw_qvalue_output.tsv", sep="\t", index=False)

def seek_pattern(pattern, uniq_save, peptide, all_ids, spectrum, text_font):
    match = re.findall(rf".*{pattern}.*", text_font, re.M|re.I)
    for result in match:
        print(result)
        result = ','.join(str(result).split("\t")[5:])
        # if matches with more than one protein, this is not an unique spectrum
        if len(result) == 0:
            # Save unique peptides in a dictionary
            # peptide is key, protein_id is value
            uniq_save.append(f'{peptide}\t{all_ids}\t{spectrum}\t')

# ---------- Parse information ----------------

# Primeira coluna
# De tras pra frente > num_charge_spectrum
# [...]/Proteome/10_TDGW20_GPI_Tripo_Tr1_Run3_1695_2_1

# open file that will be used to extract peptides
percolator_out = open(f"{args.basename}_qvalue_filtered.tsv", "r").readlines()
# open text that will be used to search peptides
percolator_search = open(f"{args.basename}_qvalue_filtered.tsv", "r").read().split("\n")
del percolator_out[0]

percolator_search = open("testes/percolator_output.tsv", "r").read().split("\n")
del percolator_search[0]
unique_peptide = []
unique_spectrum = []
# regex remover numeros estranhos no peptideo >> \[[0-9].*\]
for line in percolator_out:
    # delete line that will search (prevent false positive)
    del percolator_search[0]
    # rejoin list as text >> enable to search pattern
    percolator_search = '\n'.join(percolator_search)
    line_split = line.rstrip("\n").split("\t")
    # Protein Id > can have more than one
    all_ids = ",".join(line_split[5:])
    score = line_split[1]
    # tirar marcacoes ou procurar marcadoo?
    # regex remover numeros estranhos no peptideo >> \[[0-9].*\]
    peptide = line_split[4]
    filecontent = line_split[0].split("/")[-1]
    idname, spectrum, charge, num = filecontent.rsplit("_", 3)
    # If there are more than one protein with this peptide, this one is not unique
    if len(all_ids) == 1:
        seek_pattern(peptide, unique_peptide, peptide, all_ids, spectrum, percolator_search)
    # seek for unique spectrum
    seek_pattern(f"_{spectrum}_", unique_spectrum, peptide, all_ids, spectrum, percolator_search)
    percolator_search = percolator_search.split("\n")
    break

with open(f"{args.basename}_parsed.tsv", "r") as output:
    output.write("ProteinID\t#peptídeo único\t# de peptídeos\t# espectro unicos\t# de espectros")