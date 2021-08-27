#!/usr/bin/python3

######################################################
###   PARSE ANNOTATED FASTA FILE AND APPEND INFO   ###
###             TO NON-ANNOTATED FASTA             ###
######################################################

# USAGE: python3 script.py --annotated-fasta annotated_seq.fasta \
#                          --non-annotated-fasta non_annotated_seq.fasta \
#                         [--delimiter character]
#                         [--position integer]

# --- IMPORT MODULES -----------------------------------------------------------

import pathlib
import argparse
from Bio import SeqIO

# --- PARSER ARGUMENTS ---------------------------------------------------------

parser = argparse.ArgumentParser(
    add_help=False,
    description='''
Dad, can you make me a sandwich?

Script to parse headers from annotated fasta files.
    ''',
    epilog="""Poof, you're a sandwich.""",
    formatter_class=argparse.RawTextHelpFormatter
)

required_args = parser.add_argument_group("required arguments")
optional_args = parser.add_argument_group("custom arguments")

# mandatory arguments
required_args.add_argument(
    "-noanno", "--non-annotated-fasta", dest="noanno",
    metavar="[non-annotated file]",
    help="non-redundant fasta file",
    required=True
)

required_args.add_argument(
    "-anno", "--annotated-fasta", dest="anno",
    metavar="[annotated file]",
    help="non-redundant annotated fasta file",
    required=True
)

# optional arguments
#   no argument: uses default
optional_args.add_argument(
    "-d", "--delimiter", dest="delim",
    metavar="[annotated fasta header delimiter] (default: space delimited)",
    default=" ",
    help="delimiter for header in annotated fasta file"
)

optional_args.add_argument(
    "-pos", "--position", dest="pos",
    metavar="[id tag position] (default: [0])",
    type=int,
    default="0",
    help="id tag location after split (based on delimiter) starting from [0]"
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

# --- PREPARING SOME VARIABLES -------------------------------------------------

noanno_file = pathlib.Path(args.noanno).absolute()
noanno_filename = noanno_file.stem
anno_file = pathlib.Path(args.anno).absolute()
anno_filename = anno_file.stem

# --- FUNCTIONS ----------------------------------------------------------------

#def fasta_fetcher(input_fasta, id_list, fetcher_output):
#    wanted = sorted(set(id_list))
#    records = (seq for seq in SeqIO.parse(input_fasta, "fasta") for r in wanted if r in seq.id)
#    count = SeqIO.write(records, fetcher_output, "fasta")
#    if count < len(wanted):
#        print("IDs not found in input FASTA file")

def fasta_fetcher(input_fasta, id_list, fetcher_output):
    """find sequences in fasta file based on list containing IDs"""
    for item in id_list:
        if args.delim in item:
            item = item.split(args.delim)[args.pos]
        else:
            pass
    wanted = set(id_list)
    records = (r for r in SeqIO.parse(input_fasta, "fasta") if r.id in wanted)
    count = SeqIO.write(records, fetcher_output, "fasta")
    if count < len(wanted):
        print("IDs not found in input FASTA file")

# --- CREATE ANNOTATED FASTA ---------------------------------------------------

# file into array of lines
anno_ids = [line.strip() for line in open(args.anno) if ">" in line]
anno_ids = [sub.replace(">", "") for sub in anno_ids]
anno_tag = [tag.split(args.delim)[args.pos] for tag in anno_ids]

filtered_fasta = noanno_filename + "_filtered.fasta"
fasta_fetcher(args.noanno, anno_tag, filtered_fasta)

# annotate separated annotated products
corrected_fasta = str(noanno_filename + "_annotated_seqs.fasta")
with open(filtered_fasta) as original, open(corrected_fasta, "w") as corrected:
    records = SeqIO.parse(original, "fasta")
    for record in records:
        for item in anno_ids:
            if record.id in item:
                record.id = item
                record.description = ""
        SeqIO.write(record, corrected, "fasta")

# remove temporary file
pathlib.Path(noanno_filename + "_filtered.fasta").unlink()
