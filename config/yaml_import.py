#!/usr/bin/python3

##############################################
###   LOAD CONFIGURATION FROM .yaml FILE   ###
##############################################

# A AVENTURA VAI COMEÇAR -------------------------------------------------------

import yaml

# IMPORT CONFIG FILE -----------------------------------------------------------

with open("AnnotaTest.yaml", "r") as stream:
    try:
        config = yaml.load(stream, Loader=yaml.SafeLoader)
    except yaml.YAMLError as exc:
        print(exc)


# get element from list
#print(threads)

databases = config["databases"]
print(databases.get('swissprot-db'))

spdb_path = databases.get("secondary-db")
if str(databases.get("secondary-format")).lower() == 'eupathdb':
    flag_spdb = "-spdb"
elif str(databases.get("secondary-format")).lower() == 'trembldb':
    flag_spdb = "-trbl"
else:
    flag_spdb = "-nr"

python_exe = config['pipeline']['python']
AnnotaPipeline = config['pipeline']
databases = config['databases']
AnnotaBasename = AnnotaPipeline['basename']
keyword_list = f"'{AnnotaPipeline['keywords']}'"   # join words with '
augustus_main = config['augustus']
augustus_optional = config['augustus-optional']
seq_cleaner = config['seq-cleaner']
interpro = config['interproscan']
hmmscan = config['hmmer']
blast = config['local-aligner']
rpsblast = config['rpsblast']
kallisto = config['kallisto']
proteomics = config['proteomics']




print(f'{str(python_exe)} {str("blastp_parser.py")} -s ' \
    f' -sp ' \
    f'{str(databases.get("swissprot-db"))} -basename ' \
    f'{str(AnnotaBasename)} {str(flag_spdb)} {str(spdb_path)} -id' \
    f'{str(blast.get("identity"))} -pos' \
    f'{str(blast.get("positivity"))} -cov' \
    f'{str(blast.get("coverage"))} -kw' \
    f'{str(keyword_list)} -t' \
    f'{str(AnnotaPipeline.get("threads"))} -hsps' \
    f'{str(blast.get("hsp-max"))} -evalue {str(blast.get("evalue"))}')