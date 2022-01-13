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