#!/usr/bin/python3

# --- A AVENTURA VAI COMEÇAR ---------------------------------------------------

import argparse
import sys
from pathlib import Path
import warnings

# PARSER ARGUMENTS -------------------------------------------------------------

parser = argparse.ArgumentParser(
    add_help=False,
    description='''
Script to parse Percolator Outputs
INFO: DECOY results will also be removed
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

parcial_output = open(f"{args.basename}_qvalue_filtered.tsv", "w")
parcial_output.write(str(percolator_file[0]))
del percolator_file[0]
for line in percolator_file:
    line = line.split("\t")
    id = line[0]
    # escape decoy results
    if "DECOY_" in id:
        pass
    elif float(line[2]) <= args.qvalue:
        parcial_output.write(str('\t'.join(line)))

parcial_output.close()

# ---------- Parse information ----------------

# open file that will be used to extract peptides
percolator_out = open(f"{args.basename}_qvalue_filtered.tsv", "r").readlines()
# Check if there are valid queries for qvalue used as cutoff
if len(percolator_out) < 2:
    warnings.warn(f"INFO: file {args.basename}_qvalue_filtered.tsv return no valid queries for qvalue", stacklevel=2)
    sys.exit(2)
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
        # If list contains only one element, get this guy
        only_id = all_ids[0]
        if only_id not in dict_write.keys():
            dict_write[only_id] = [f"{peptide}\t{spectrum}\n"]
        else:
            dict_write[only_id].append([f"{peptide}\t{spectrum}\n"])


with open(f"{args.basename}_parsed.tsv", "w") as output:
    output.write("ProteinID\tPeptide\tSpectrum\n")
    for key in sorted(dict_write.keys()) :
        info = dict_write.get(key)
        for peptide in info:
            output.write(f"{key}\t{''.join(peptide)}")
