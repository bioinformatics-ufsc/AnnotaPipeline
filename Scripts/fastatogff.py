#!/usr/bin/python3

import argparse

# --- PARSER ARGUMENTS ---------------------------------------------------------

parser = argparse.ArgumentParser(
        add_help=False,
        description='''
        AnnotaPipeline script to parse GFF

        Input GFF file and All annotated Products to join annotation into GFF

        NOTE: If you make a subset from All_Annotated_Products.txt and input in this script,
        the output will have only this subset.
        ''',

        epilog="""Hear the whispers of your hope, the answer wasn't told""",
        formatter_class=argparse.RawTextHelpFormatter
)

requiredNamed = parser.add_argument_group('required arguments')
optionalNamed = parser.add_argument_group('optional arguments')

# mandatory arguments
#   type (default): string
requiredNamed.add_argument(
        '-all', '--allah', dest='anot',
        metavar='[All_Annotated_Products.txt]',
        help='One of AnnotaPipeline Outputs',
        required=True
)

requiredNamed.add_argument(
        '-gff', '--gff', dest='gff',
        metavar='[GFF file]',
        help='GFF file with all predictions',
        required=True
)

requiredNamed.add_argument(
        '-b', '--basename', dest='basename',
        metavar='[Jonas]',
        help='Just a name to your child',
        required=True
)

# custom [--help] argument
optionalNamed.add_argument(
        '-h', '-help', '--help',
        action='help',
        default=argparse.SUPPRESS,
        help='What kind of sorcery is this?'
)

# arguments saved here
args = parser.parse_args()

annotations = open(args.anot, "r")
anots = []
ids = []
# Get annotations and ids
while True:
        line = annotations.readline()
        array_anot = line.split("\t")
        try:
                anots.append(array_anot[1])
                ids.append(array_anot[0])
        except Exception as warn:
                pass
        if not line:
                break

# Build dictionary with id and anotation
dict_anot = dict(zip(ids, anots))
del ids, anots

gff = open(str(args.gff), "r").read().split("# start gene ")
header_gff = gff[0]
del gff[0]

outfile = open(f"{str(args.basename)}_Annotated_GFF.gff", "w")
outfile.write(str(header_gff))
infos = []

for gene in gff:
        gene = gene.split("\n")
        gene_stats = gene[1:]
        infos.append(f"# start gene {str(gene[0])}")
        id_gene = None
        count = 0
        count_transcript = None
        transcript_line = None
        for stat in gene_stats:
                count += 1
                if "\ttranscript\t" not in stat:
                        infos.append(stat)
                else:
                        stat_split = stat.split("\t")
                        #id_gene = stat_split[-1]
                        # use regex to remove ID= and ;PARENT from annotation ID in gff AUGUSTUS output
                        id_gene =  re.sub(r';.*', '', stat_split[-1])
                        id_gene =  re.sub(r'.*=', '', id_gene)
                        count_transcript = count
                        transcript_line = "\t".join(stat_split[0:-1])
                if stat is gene_stats[-1]:
                        if id_gene in dict_anot.keys():
                                try:
                                        outfile.write("\n".join(infos[0:count_transcript]) + "\n" +
                                                                  transcript_line + "\t" + dict_anot.get(id_gene)
                                                                  + "\n".join(infos[count_transcript:]))
                                except NameError:
                                        raise Exception("Inconsistent ID between gff and annotations")
                        infos = []
                        count_transcript = None
                        transcript_line = None
