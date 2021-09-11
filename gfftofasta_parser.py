#!/usr/bin/env python3

import argparse
import warnings

parser = argparse.ArgumentParser(
    add_help=False,  # removes original [--help]
    description='''Script to write gff info into fasta sequences
Squence's headers will be like this:
><sequence_id> | Organism: <given organism> | Location: <scaffold/chromossome> | Start: <info> | End: <info> | Strand: <+/-> | Description: <info>
Example:
>g1.t1 | Organism: Homo sapiens | Location: scaffold1 | Start: 12768 | End: 14450 | Strand: + | Description: Glucose-6-phosphate 1
-dehydrogenase 5, cytoplasmic		 ''',
    epilog="""We stand before the dawn of a new world.""", formatter_class=argparse.RawTextHelpFormatter
)

requiredNamed = parser.add_argument_group('required arguments')
optionalNamed = parser.add_argument_group('optional arguments')

# mandatory arguments
#   type (default): string
requiredNamed.add_argument(
    '-gff', dest='gff',
    metavar='[augustus_prediction.gff]',
    help='Augustus prediction file',
    required=True
)

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

# Info from transcript products will be stored here
trans = {}
# Info from annotation will be stored here
ant = {}

# ============================= Extract info from gff =========================
gff = open(str(args.gff), "r").read().split("# start gene ")
# Delete header with information of model and other blablabla
del gff[0]

for gene in gff:
    # Get second line - transcript atributes
    gene = gene.split("\n")
    gene_features = gene[2].split()
    # Save all info into a dictionary
    trans[gene_features[-1]] = {"scaff": str(gene_features[0]), "start": str(gene_features[3]),
                                "end": str(gene_features[4]), "strand": str(gene_features[6])}

# ============================= Extract annotations ==========================

annot = open(str(args.annot)).read().splitlines()
for line in annot:
    line = line.split("\t")
    ant[line[0]] = line[1]

fasta = open(str(args.seqs), "r").read().split(">")
del fasta[0]  # remove empty value

# ============== Write fasta file with new header ==========================

# Remove quotes from organism 
organism = str(args.org).replace("\"", "").replace("\'", "")

# annotafile = open("AnnotaPipeline_" + str(args.basename) + ".fasta", "w")
with open(f"AnnotaPipeline_{str(args.basename)}_proteins.fasta", "w") as annotafile:
    ids_warn = []  # store ids with no annotation
    for sequence in fasta:
        sequence = sequence.split("\n")
        id = sequence[0]
        try:
            annotafile.write(
                f">{str(id)} | "
                f"Organism: {organism} | "
                f"Location: {str(trans[id]['scaff'])} | "
                f"Start: {str(trans[id]['start'])} | "
                f"End: {str(trans[id]['end'])} | "
                f"Strand: {str(trans[id]['strand'])} | "
                f"Description: {str(ant[id])}\n"
            )
            annotafile.write(str("\n".join(sequence[1:])))
        except Exception as warn:
            # warn ignored, user will be informed all at once
            ids_warn.append(id)
            pass

# annotafile.close()

# =============== Warning message, if something goes wrong =====================

if len(ids_warn) > 0:
    warnings.warn("WARNING: Not all sequences from fasta file were in annotation file", stacklevel=2)
    warn_seq = open(f"{str(args.basename)}_ids_with_no_annotations.txt", "w")
    for a in ids_warn:
        warn_seq.write(str(a) + "\n")
    warn_seq.close()
    warnings.warn("INFO: IDs stored in ids_with_no_annotations.txt", stacklevel=2)
