#!/usr/bin/python3

# A AVENTURA VAI COMEÇAR -------------------------------------------------------

import yaml

# IMPORT CONFIG FILE -----------------------------------------------------------

with open("AnnotaPipeline.yaml", "r") as stream:
    try:
        yaml_dict = yaml.load(stream, Loader=yaml.SafeLoader)
    except yaml.YAMLError as exc:
        print(exc)

print(yaml_dict)

# get element from list
threads = yaml_dict["pipeline"]["threads"]
print(threads)