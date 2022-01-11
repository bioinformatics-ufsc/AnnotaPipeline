#!/usr/bin/env python3

import argparse
import warnings

parser = argparse.ArgumentParser(
	add_help=False,  # removes original [--help]
	description='''Script to join sequence and annotations in a simple way (alternative to gff_to fasta, if you don't
have a gff file)
Squence's headers will be like this:
><sequence_id> | Organism: <given organism> | Description: <info>
Example:
>g1.t1 | Organism: Homo sapiens | Description: Glucose-6-phosphate 1
-dehydrogenase 5, cytoplasmic		 ''',
	epilog="""We stand before the dawn of a new world.""", formatter_class=argparse.RawTextHelpFormatter
)

requiredNamed = parser.add_argument_group('required arguments')
optionalNamed = parser.add_argument_group('optional arguments')

# mandatory arguments
#   type (default): string

requiredNamed.add_argument(
	'-annot', '--annotations', dest='annot',
	metavar='[annotations.txt]',
	help='File with id and annotation separated by \t',
	required=True
)

requiredNamed.add_argument(
	'-b', '--basename', dest='basename',
	metavar='[Jonas]',
	help='It\'s a boy, and will be called Jonas',
	required=True
)

requiredNamed.add_argument(
	'-faf', '--fasta', dest='seqs',
	metavar='[fasta.fasta]',
	help='Fasta file with protein sequences, extracted from augustus gff',
	required=True
)

requiredNamed.add_argument(
	'-org', '--organism', dest='org',
	metavar='[Homo sapiens]',
	help='Organism name. PS: if you want to pass genus and species give it with " - ex: "Homo Sapiens"',
	required=True
)

# custom [--help] argument
optionalNamed.add_argument(
	'-h', '-help', '--help',
	action='help',
	default=argparse.SUPPRESS,  # hidden argument
	help='It\'s going to be legen - wait for it - dary!'
)

# --------------------------------------------------------

# arguments saved here
args = parser.parse_args()

# Info from annotation will be stored here
ant = {}

# ============================= Extract annotations ==========================

annot = open(str(args.annot)).read().splitlines()
for line in annot:
	line = line.split("\t")
	ant[line[0]] = line[1]

fasta = open(str(args.seqs), "r").read().split(">")
del fasta[0]  # remove empty value

# ============== Write fasta file with new header ==========================

out = open(f"AnnotaPipeline_{str(args.basename)}_proteins.fasta", "w")
ids_warn = []  # store ids with no annotation
for sequence in fasta:
	sequence = sequence.split("\n")
	id = sequence[0]
	try:
		out.write(f">{str(id)} | Organism: {str(args.org)} | Description: {str(ant[id])}\n")
		out.write(str("\n".join(sequence[1:])))
	except Exception as warn:
		ids_warn.append(id)
		pass

out.close()

# =============== Warning message, if something goes wrong =====================

if len(ids_warn) > 0:
	warnings.warn("WARNING: Not all sequences from fasta file were in annotation file", stacklevel=2)
	warn_seq = open(f"{str(args.basename)}_ids_with_no_annotations.txt", "w")
	for a in ids_warn:
		warn_seq.write(str(a) + "\n")
	warn_seq.close()
	warnings.warn("INFO: IDs stored in ids_with_no_annotations.txt", stacklevel=2)