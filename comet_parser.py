#!/usr/bin/python3

# --- A AVENTURA VAI COMEÇAR ---------------------------------------------------

import argparse
import pandas as pd
from pathlib import Path

# PARSER ARGUMENTS -------------------------------------------------------------

parser = argparse.ArgumentParser(
    add_help=False,
    description='''
Script to parse Comet output (tsv file)

Output given must have the following columns:

scan    num     charge  exp_neutral_mass        calc_neutral_mass       e-value xcorr   
delta_cn        sp_score        ions_matched      ions_total      plain_peptide   modified_peptide
prev_aa next_aa protein protein_count   modifications
    ''',
    epilog="""Shadows will fade... someday...""",
    formatter_class=argparse.RawTextHelpFormatter
)

required_args = parser.add_argument_group("required arguments")
optional_args = parser.add_argument_group("custom arguments")

# mandatory arguments
required_args.add_argument(
    "-p", "--proteom-file", dest="proteom",
    type=argparse.FileType('r'), 
    metavar="[Comet output]",
    help="Peptide quantification file (Comet output)",
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
optional_args.add_argument(
    "-ch", "-charge", dest='charge',
    metavar = "[2]",
    type=int, 
    help="Parse output also by column 'charge'. Higher or equal than [value]",
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

# --- PARSE DATAFRAMES ---------------------------------------------------------

df = pd.read_csv(args.proteom, sep="\t", skiprows=1, index_col=False)
with pd.option_context('mode.chained_assignment', None):
    # Delete CometVersion and header lines inside file
    df_no_comet_version = df[df["scan"].str.contains("CometVersion|scan", regex=True, na=False)==False]
    df_no_comet_version['charge'] = df_no_comet_version['charge'].apply(pd.to_numeric)
    # Create decoy_true and decoy_false datasets
    decoy_false = df_no_comet_version[df_no_comet_version["protein"].str.contains("DECOY_", regex=False, na=False)==False]
    decoy_true = df_no_comet_version[df_no_comet_version["protein"].str.contains('DECOY_', regex=False, na=False)==True]

# If needed, filter by charge
if args.charge is not None:
    with pd.option_context('mode.chained_assignment', None):
        decoy_false = decoy_false.loc[decoy_false['charge'] >= int(args.charge)]
        decoy_true = decoy_true.loc[decoy_true['charge'] >= int(args.charge)]


decoy_false.to_csv(f"{args.basename}_peptide_identificated_DECOY_FALSE.txt", sep="\t", index=False)
decoy_true.to_csv(f"{args.basename}_peptide_identificated_DECOY_TRUE.txt", sep="\t", index=False)