from setuptools import setup, find_packages

requirements = [
                'biopython>=1.73',
                'pandas>=0.24.1',
                'pyyaml>=6.0'
                #'shutil',
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
    data_files=[('config', ['config/AnnotaPipeline.config'])],
    scripts=[
                'Scripts/AnnotaPipeline.py',
                'Scripts/blastp_parser.py',
                'Scripts/fasta_simple.py',
                'Scripts/fastatogff.py',
                'Scripts/funcannotation_parser.py',
                'Scripts/gfftofasta_parser.py',
                'Scripts/info_parser.py',
                'Scripts/kallisto_parser.py',
                'Scripts/percolator_parser.py',
                'Scripts/summary_parser.py'
             ],
    include_package_data=True,
    python_requires=">=3.6.9",
    url='https://github.com/GuiMaia/AnnotaPipeline.git',
    download_url='https://github.com/GuiMaia/AnnotaPipeline/archive/refs/heads/conda_env.zip',
# USA esse cara pra fz o classifier
# https://pypi.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: Apache Software License",
        "Intended Audience :: Science/Research",
        "Natural Language :: Portuguese (Brazilian)"
        "Operating System :: Unix",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9"
    ],
)
