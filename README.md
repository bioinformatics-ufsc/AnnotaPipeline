[![Quote](https://img.shields.io/badge/Quote-May%20the%20Dark%20shine%20your%20way%20(Darkdiver%20Grandahl)-blueviolet)](https://g.co/kgs/bMJLfj)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

<font size=16>**AnnotaPipeline**</font>

- [**About AnnotaPipeline**](#about-annotapipeline)
- [**Individual analysis**](#individual-analysis)
  - [**1. Gene prediction**](#1-gene-prediction)
  - [**2. Similarity analysis**](#2-similarity-analysis)
  - [**3. Functional annotation**](#3-functional-annotation)
  - [**4. Transcript quantification**](#4-transcript-quantification)
  - [**5. Peptide idenfication**](#5-peptide-idenfication)
- [**How to setup**](#how-to-setup)
  - [**Download databases**](#download-databases)
  - [**Using conda**](#using-conda)
    - [**Setup AUGUSTUS species for personalized predictions**](#setup-augustus-species-for-personalized-predictions)
  - [**Manual install**](#manual-install)
  - [**Setup `AnnotaPipeline.yaml`**](#setup-annotapipelineyaml)
- [**How to run**](#how-to-run)
- [**Output**](#output)


# **About AnnotaPipeline**

Integrated tool to annotate hypothetical proteins developed by the Laboratório de Bioinformática at Universidade Federal de Santa Catarina (Brazil).

If you have questions, suggestions or difficulties regarding the pipeline, please do not hesitate to contact our team here on GitHub or by email: <labioinfo.genome@gmail.com>.


# **Individual analysis**

## **1. Gene prediction**

## **2. Similarity analysis**

## **3. Functional annotation**

## **4. Transcript quantification**

## **5. Peptide idenfication**


# **How to setup**

AnnotaPipeline requires the following software to run properly:
  * BLAST+ and RPS BLAST (available at <https://ftp.ncbi.nih.gov/blast/executables/blast+/LATEST>)
  * InterProScan (available at <https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html>)
  * HMMER (available at <http://hmmer.org/download.html>)

You will also need to install [AUGUSTUS](https://doi.org/10.1093/bioinformatics/btr010) (available at <https://github.com/Gaius-Augustus/Augustus>) if you want to run this pipeline starting with gene/protein prediction.

Before executing, please modify the necessary fields in the configuration file (`AnnotaPipeline.yaml`).

## **Download databases**

**Required:**
  * SwissProt (avaliable at <https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz>)
  * Pfam (avaliable at <http://pfam.xfam.org>)
  * CDD (avaliable at <https://www.ncbi.nlm.nih.gov/cdd>)

**Choose one secondary database:**
  * TrEMBL (avaliable at <https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz>)
  * EuPathDB (avaliable at <https://veupathdb.org/veupathdb/app>):
    * AmoebaDB
    * CryptoDB
    * FungiDB
    * GiardiaDB
    * HostDB
    * MicrosporidiaDB
    * PiroplasmaDB
    * PlasmoDB
    * ToxoDB    (tested)
    * TrichDB
    * TriTrypDB (tested)
  * NR Database | NCBI (available at <https://ftp.ncbi.nlm.nih.gov/blast/db>)
  > **TIP:** You can use a subset of NR Database

## **Using conda**

1. Download `conda_environment.yaml` file

2. Create envinroment
```bash
    conda env create -n <desired_name> -f conda_environment.yaml
```

3. Activate environment
```bash
    conda activate <desired_name>
```

### **Setup AUGUSTUS species for personalized predictions**

1. Locate AnnotaPipeline environment home
```bash
    echo $CONDA_PREFIX
```

2. Go to `$CONDA_PREFIX/config/species`
```bash
    cd $CONDA_PREFIX/config/species
```

3. Add custom species folder (trained results)

## **Manual install**

1. Clone repository
```bash
    git clone https://github.com/GuiMaia/AnnotaPipeline.git
```

2. Run `setup.py` (Scripts will be avaliable at `$PATH`)

2. Install required softwares:
  * BLAST+ and RPS BLAST (available at <https://ftp.ncbi.nih.gov/blast/executables/blast+/LATEST>)
  * InterProScan (available at <https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html>)
  * HMMER (available at <http://hmmer.org/download.html>)

3. Optional softwares
  * Kallisto (available at <https://pachterlab.github.io/kallisto/download.html>)
  * Comet MS/MS (available at <https://github.com/UWPR/Comet/releases/latest>)
    * Requires Percolator (available at <https://github.com/percolator/percolator>)

4. Download databases
5. Configure `AnnotaPipeline.yaml`

## **Setup `AnnotaPipeline.yaml`**
> **TIP:** If you already have InterProScan locally installed and configured, use it instead (as an alternative to conda installation &ndash; interpro vanilla)

We recommend using our example configuration file as a guide (`config_example.yaml`).


# **How to run**

To initiate the pipeline, you can run `AnnotaPipeline.py` with three different options:

```bash
python3 AnnotaPipeline.py -p protein_sequences.fasta
```

This is the most simple version of AnnotaPipeline execution.

The annotation process will begin with the submitted `protein_sequences.fasta` and an simplified version of `BASENAME_Annotated_GFF.gff` file at the end.

```bash
python3 AnnotaPipeline.py -s genomic_data.fasta
```

This is the most complete version of AnnotaPipeline execution.


It will execute gene and protein prediction on the `genomic_data.fasta` utilizing AUGUSTUS, and then utilize the predicted proteins to initiate the annotation process.

Because of this first step, it is important that you have trained AUGUSTUS to your particular species before executing AnnotaPipeline.

```bash
python3 AnnotaPipeline.py -p protein_sequences.fasta -gff gff_file.gff
```

You can execute the AnnotaPipeline utilizing this command line, in case you already have a protein sequence file and a GFF3 file from previous AUGUSTUS predictions.

The annotation process is the same as the first option, the difference being you will also have an annotated GFF3 file as output.

Please notice that the submited ```gff_file.gff``` needs to be in GFF3 format.


# **Output**
Depending on which option you decided to execute the AnnotaPipeline, it will output three files (among many others in their respective folders):

* `All_Annotated_Products.txt` - which contains all unique sequence identifiers and their respective annotations (with functional annotations, when present).
* ` Annota_BASENAME.fasta` - which contains all the sequences and their annotations (with functional annotations, when present) in FASTA format.
* `BASENAME_Annotated_GFF.gff` - which contains all the sequences and their annotations (with functional annotations, when present) in GFF3 format.
* `AnnotaPipeline_BASENAME_transcripts.fasta` - which contains nucleotide sequence for predicted proteins, with same features present at protein file.
* `AnnotaPipeline_BASENAME_Summary.tsv` - which summarize hits found for each protein, in similarity, functional, transcriptomics (if used) and proteomics analysis (if used).

Raw outputs are listed inside output folder:

* `1_GenePrediction_BASENAME`           - AUGUSTUS files
* `2_SimilarityAnalysis_BASENAME`       - BLAST analysis
* `3_FunctionalAnnotation_BASENAME`     - InterProscan/HMM/RPSBLAST analysis
* `4_TranscriptQuantification_BASENAME` - Kallisto analysis
* `5_PeptideIdentification_BASENAME`    -

The output folders ysisand files will be located in the same folder you executed the pipeline in.