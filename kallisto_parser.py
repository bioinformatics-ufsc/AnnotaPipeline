#!/usr/bin/python3

# A AVENTURA VAI COMEÇAR -------------------------------------------------------

from pathlib import Path
import argparse
import pandas as pd

# PARSER ARGUMENTS -------------------------------------------------------------

parser = argparse.ArgumentParser(
    add_help=False,
    description='''
Dad, can you make me a sandwich?
    ''',
    epilog="""Poof, you're a sandwich.""",
    formatter_class=argparse.RawTextHelpFormatter
)

required_args = parser.add_argument_group("required arguments")
optional_args = parser.add_argument_group("custom arguments")

# mandatory arguments
required_args.add_argument(
    "-ktfile", "--kallisto-file", dest="ktfile",
    metavar="[Kallisto output]",
    help="Abundance TSV file (Kallisto output)",
    required=True
)

required_args.add_argument(
    "-tpm", "--tpm-threshold", dest="tpm",
    metavar="[TPM threshold]",
    help="TPM threshold value to filter Kallisto output",
    required=True
)

# custom [--help] argument
optional_args.add_argument(
    "-h", "-help", "--help",
    action="help",
    default=argparse.SUPPRESS,
    help="What kind of sorcery is this?"
)

# arguments saved here
args = parser.parse_args()

# SET FILE LOCATION ------------------------------------------------------------

files_dir = Path("/mnt/bfe2b7ea-8d6d-4d3e-ae29-96991c6abfaf/mestrado/AnnotaPipeline")

# FILES ------------------------------------------------------------------------

kallisto_file = Path(files_dir / "abundance.tsv")

# OPEN DATAFRAMES --------------------------------------------------------------

df_kallisto = pd.read_csv(kallisto_file, sep="\t")

# PARSE DATAFRAMES -------------------------------------------------------------

tpm_avg = df_kallisto["tpm"].mean()
tpm_median = df_kallisto["tpm"].median()

df_avg = df_kallisto.loc[df_kallisto["tpm"] >= tpm_avg]
df_median = df_kallisto.loc[df_kallisto["tpm"] >= tpm_median]