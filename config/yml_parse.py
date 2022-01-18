import yaml

# with open("AnnotaPipeline.yaml", "r") as stream:
#     try:
#         print(yaml.safe_load(stream))
#     except yaml.YAMLError as exc:
#         print(exc)


a_yaml_file = open("AnnotaPipeline.yaml")
parsed_yaml_file = yaml.load(a_yaml_file, Loader=yaml.FullLoader)

# Get list
pipeline = parsed_yaml_file.get("pipeline")
pip = parsed_yaml_file["pipeline"]
print(pipeline)
print(pip)

# Get element from list
threads = parsed_yaml_file["pipeline"].get('threads')
print(threads)

