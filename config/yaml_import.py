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

print(config)

# get element from list
threads = config["pipeline"]["threads"]
print(threads)