from setuptools import setup, find_packages

# setar as versoes sertinho (>=x ou ==x)
requirements = [
                'bioopython>=1.73',
                #'shutil',
                'pandas>=0.24.1',
                #'argparse',
                #'configparser',
                #'logging',
                #'os',
                #'pathlib',
                #'subprocess',
                #'sys',
                #'re',
                #'warnings'
                ]

# default python libraries
# https://docs.python.org/3/library/

with open('README.md') as rm:
    long_description = rm.read()

setup(
    name="AnnotaPipeline",
    packages=find_packages(),
    version='1.0',
    description='Pipeline to predict and annotate eukaryotic genes',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=requirements,
    author='Guilherme Augusto Maia, Eric Kazuo Kawagoe, Vilmar Benetti Filho, \
            Tatiany Aparecida Soratto, Renato Moreira Simões, Glauber Wagner',
    author_email='labinfo.ufsc@gmail.com',
    zip_safe=False,
    data_files=[('config', ['AnnotaPipeline.config'])],
    scripts=['AnnotaPipeline.py',
             'blastp_parser.py',
             'fasta_simple.py',
             'fastatogff.py',
             'funcannotation_parser.py',
             'gfftofasta_parser.py',
             'info_parser.py',
             'kallisto_parser.py',
             'percolator_parser.py',
             'summary_parser.py'
             ], # adicionar cada um dos scripts
    include_package_data=True,
    python_requires=">=3.6.9", # versao minima pra ter pathlib
    url='https://github.com/GuiMaia/AnnotaPipeline/tree/v1.0.git',
    download_url='https://github.com/GuiMaia/AnnotaPipeline/archive/refs/heads/v1.0.zip',
    classifiers=[
    ],
)

# Pangolin exemple
# from setuptools import setup, find_packages
# import glob
# import os
# import pkg_resources

# from pangolin import __version__, _program

# setup(name='pangolin',
#       version=__version__,
#       packages=find_packages(),
#       scripts=['pangolin/scripts/pangolearn.smk'
#                 ],
#       package_data={"pangolin":["data/*"]},
#       install_requires=[
#             "biopython>=1.70",
#             'pandas>=1.0.1',
#             "wheel>=0.34",
#             'joblib>=0.11',
#             'scikit-learn>=0.23.1',
#             "PuLP>=2"
#         ],
#       description='phylogenetic assignment of named global outbreak lineages',
#       url='https://github.com/cov-lineages/pangolin',
#       author='Aine OToole & Emily Scher',
#       author_email='aine.otoole@ed.ac.uk',
#       entry_points="""
#       [console_scripts]
#       {program} = pangolin.command:main
#       """.format(program = _program),
#       include_package_data=True,
#       keywords=[],
#       zip_safe=False)