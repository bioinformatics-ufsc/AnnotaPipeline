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
keyword_list = AnnotaPipeline['keywords']   # join words with '
augustus_main = config['augustus']
augustus_optional = config['augustus-optional']
seq_cleaner = config['seq-cleaner']
interpro = config['interproscan']
hmmscan = config['hmmer']
blast = config['local-aligner']
rpsblast = config['rpsblast']
kallisto = config['kallisto']
proteomics = config['proteomics']

kw = str(f'{",".join(keyword_list)}')
print(kw)
k2 = kw.split(",")
for word in k2:
    print(word)
