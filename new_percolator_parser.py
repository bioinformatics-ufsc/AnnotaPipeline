#!/usr/bin/python3

# --- A AVENTURA VAI COMEÇAR ---------------------------------------------------

import argparse
import sys
from pathlib import Path
import pandas as pd
import warnings
from collections import defaultdict

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



percolator_file = open(str(args.proteom), "r").readlines()
# If needed, filter by charge
# PANDAS N FUNCIONAAAAAAAAAAAAA A
parcial_output = open(f"{args.basename}_qvalue_filtered.tsv", "w")
parcial_output.write(str(percolator_file[0]))
del percolator_file[0]
for line in percolator_file:
    line = line.split("\t")
    if float(line[2]) <= args.qvalue:
        parcial_output.write(str('\t'.join(line)))

parcial_output.close()

# ---------- Parse information ----------------

# Primeira coluna
# De tras pra frente > num_charge_spectrum
# [...]/Proteome/10_TDGW20_GPI_Tripo_Tr1_Run3_1695_2_1

# open file that will be used to extract peptides
percolator_out = open(f"{args.basename}_qvalue_filtered.tsv", "r").readlines()
# Check if there are valid queries for qvalue used as cutoff
if len(percolator_out) < 2:
    warnings.warn(f"INFO: file {args.basename}_qvalue_filtered.tsv return no valid queries for qvalue", stacklevel=2)
    sys.exit(2)
# open text that will be used to search peptides
# Remove index
del percolator_out[0]

dict_write = {}
# regex remover numeros estranhos no peptideo >> \[[0-9].*\]
for line in percolator_out:
    line_split = line.rstrip("\n").split("\t")
    # Protein Id > can have more than one
    all_ids = line_split[5:]
    peptide = line_split[4]
    filecontent = line_split[0].split("/")[-1]
    idname, spectrum, charge, num = filecontent.rsplit("_", 3)
    if len(all_ids)> 1:
        for proteinID in all_ids:
            if proteinID not in dict_write.keys():
                dict_write[f'{proteinID}'] = [f"{peptide}\t{spectrum}\n"]
            else:
                dict_write[f'{proteinID}'].append([f"{peptide}\t{spectrum}\n"])
    else:
        if all_ids not in dict_write.keys():
            dict_write[all_ids] = [f"{peptide}\t{spectrum}\n"]
        else:
            dict_write[all_ids].append([f"{peptide}\t{spectrum}\n"])
    break

with open(f"{args.basename}_parsed.tsv", "w") as output:
    output.write("ProteinID\tPeptide\tSpectrum\n")
    for key in sorted(dict_write.keys()) :
        info = dict_write.get(key)
        for peptide in info:
            output.write(f"{key}\t{''.join(peptide)}")

def quantitative_proteomics(path):
    parsed_files = Path(path).glob('*_parsed.tsv')
    # Create basal 
    data = pd.DataFrame({'ProteinID':[], \
                        'Peptide':[], \
                        'Spectrum':[]}) 
    for file in parsed_files:
        df = pd.read_csv(file, sep='\t', header=0)
        data = data.append(df)
    
    # separar os duplicados
    #duplicated_pep = data[data.duplicated(['Peptide'], keep=False)]
    #unique_pep = data[data.duplicated(['Peptide'], keep=False)==False]

    #duplicated_spec = data[data.duplicated(['Peptide'], keep=False)]
    #unique_spectrum = data[data.duplicated(['Spectrum'], keep=False)==False]

    total = pd.DataFrame({'ProteinID':[], \
                    'Unique Peptide':[], \
                    'Total Peptide':[], \
                    'Unique Spectrum':[], \
                    'Total Spectrum':[]
                    }) 
    # REPETIR PARA SPECTRUM
    parcial = data.groupby(['ProteinID', 'Peptide']).size().sort_values(ascending=False).reset_index(name='Unique_Peptide')
    parcial2 = parcial.groupby(['ProteinID']).size().sort_values(ascending=False).reset_index(name='Total')
    tmp = parcial.loc[parcial['Unique_Peptide'] == 1].drop(columns=["Peptide"])
    tmp2 = tmp.groupby("ProteinID").count().reset_index()
    total["ProteinID"] = parcial2["ProteinID"]
    total.set_index("ProteinID").join(tmp2.set_index("ProteinID")).reset_index()
    # Tem q dar um sort talvez? ou usar join que garante que pega o proteinID? sla
    total["Total Peptide"] = parcial2["Total"]

    #parcial.loc[parcial['Total_Peptide'] == 1]
    #data.groupby(['ProteinID', 'Spectrum']).size().sort_values(ascending=False).reset_index(name='Total_Spectrum')

    # Fazer as contas
    #data.groupby('ProteinID').count()
    #data.groupby('Peptide').count()
    #data.groupby(['ProteinID', 'Peptide']).count()
    # THHIIISSS

    # escrever sem o index
    #data.to_csv(filename, index=False)
    print(data)
