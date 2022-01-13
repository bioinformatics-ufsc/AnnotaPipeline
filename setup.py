from setuptools import setup, find_packages

# setar as versoes sertinho (>=x ou ==x)
requirements = [
                'Bio',
                'shutil',
                'pandas',
                'argparse',
                'configparser',
                'logging',
                'os',
                'pathlib',
                'shutil',
                'subprocess',
                'sys',
                're',
                'warnings'
                ]

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
            Tatiany Aparecida Soratto, Renato Simões Glauber Wagner',
    author_email='labinfo.ufsc@gmail.com',
    zip_safe=False,
    package_data={'': ['*.config']},
    include_package_data=True,
    python_requires=">=3",
    url='https://github.com/GuiMaia/AnnotaPipeline/tree/v1.0',
    download_url='https://github.com/GuiMaia/AnnotaPipeline/archive/refs/heads/v1.0.zip',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',  # pathlib is born
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'
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