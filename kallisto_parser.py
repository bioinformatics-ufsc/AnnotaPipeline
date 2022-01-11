#!/usr/bin/python3

# A AVENTURA VAI COMEÇAR -------------------------------------------------------

from pathlib import Path
import argparse
import pandas as pd

# PARSER ARGUMENTS -------------------------------------------------------------

parser = argparse.ArgumentParser(
    add_help=False,
    description='''
Script to parse Kallisto output based on TPM threshold
    ''',
    epilog="""Default TPM threshold is 1""",
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
    "-basename", dest="basename",
    metavar="[It\'s a boy, and will be called Jonas]",
    help="basename",
    required=True
)

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument(
	"-tpmavg", "--tpm-average",  dest="tpm_avg", action="store_true",
	help='use TPM average as threshold'
)

group.add_argument(
	"-tpmmd", "--tpm-median", dest="tpm_median", action="store_true",
	help="use TPM median as threshold"
)

group.add_argument(
	"-tpmval", "--tpm-threshold", dest="tpm_value",
	metavar="-tpmth 100",
	help="use a specific number as threshold for tpm"
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

# SET FILE LOCATION ------------------------------------------------------------

kallisto_file = Path(args.ktfile).absolute()

# OPEN DATAFRAMES --------------------------------------------------------------

df_kallisto = pd.read_csv(kallisto_file, sep="\t")

# PARSE DATAFRAMES -------------------------------------------------------------

if args.tpm_avg == True:
    # remove SettingWithCopyWarning for next operation:
    with pd.option_context('mode.chained_assignment', None):
        tpm = df_kallisto["tpm"].mean()
        print(f"Average TPM: {tpm}")
        df_kt_threshold = df_kallisto.loc[df_kallisto["tpm"] >= float(tpm)]
        df_kt_threshold.drop(columns=["length", "eff_length", "est_counts"], inplace=True)
elif args.tpm_median == True:
    # remove SettingWithCopyWarning for next operation:
    with pd.option_context('mode.chained_assignment', None):
        tpm = df_kallisto["tpm"].median()
        print(f"Median TPM: {tpm}")
        df_kt_threshold = df_kallisto.loc[df_kallisto["tpm"] >= float(tpm)]
        df_kt_threshold.drop(columns=["length", "eff_length", "est_counts"], inplace=True)
else:
    # remove SettingWithCopyWarning for next operation:
    with pd.option_context('mode.chained_assignment', None):
        tpm = float(args.tpm_value)
        print(f"Chosen TPM value: {tpm}")
        df_kt_threshold = df_kallisto.loc[df_kallisto["tpm"] >= float(tpm)]
        df_kt_threshold.drop(columns=["length", "eff_length", "est_counts"], inplace=True)

df_kt_threshold.to_csv(f"{args.basename}_Transcript_Quantification.tsv", sep="\t", index=False)