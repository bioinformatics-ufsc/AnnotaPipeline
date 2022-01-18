#!/usr/bin/python3

##############################################
###   LOAD CONFIGURATION FROM .yaml FILE   ###
##############################################

# A AVENTURA VAI COMEÇAR -------------------------------------------------------

import yaml

# IMPORT CONFIG FILE -----------------------------------------------------------

with open("AnnotaPipeline.yaml", "r") as stream:
    try:
        config = yaml.load(stream, Loader=yaml.SafeLoader)
    except yaml.YAMLError as exc:
        print(exc)


# get element from list
threads = config["pipeline"]["threads"]


for section, list_section in config.items():
    if str(section) == "databases":
        if config['databases']['secondary-format'] not in ('eupathdb', 'nrdb', 'trembldb'):
            print("jonas")
        if config['databases']['secondary-format'] is None:
            print("aaaa")
    elif str(section) == "comet":
        pass