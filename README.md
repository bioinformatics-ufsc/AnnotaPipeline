[![Quote](https://img.shields.io/badge/Quote-May%20the%20Dark%20shine%20your%20way%20(Darkdiver%20Grandahl)-blueviolet)](https://g.co/kgs/bMJLfj)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

<font size=16>**AnnotaPipeline**</font>

- [**About AnnotaPipeline**](#about-annotapipeline)
- [**Analysis description**](#analysis-description)
  - [**1. Gene prediction**](#1-gene-prediction)
  - [**2. Similarity analysis**](#2-similarity-analysis)
  - [**3. Functional annotation**](#3-functional-annotation)
  - [**4. Transcript quantification (Optional)**](#4-transcript-quantification-optional)
  - [**5. Peptide idenfication (Optional)**](#5-peptide-idenfication-optional)
  - [**Last things last**](#last-things-last)
- [**How to setup**](#how-to-setup)
  - [**Download databases**](#download-databases)
  - [**Using conda**](#using-conda)
    - [**Setup AUGUSTUS species for personalized predictions**](#setup-augustus-species-for-personalized-predictions)
  - [**Manual install**](#manual-install)
  - [**Setup `AnnotaPipeline.yaml`**](#setup-annotapipelineyaml)
- [**How to run**](#how-to-run)
  - [**Protein file as input**](#protein-file-as-input)
  - [**Genome file as input**](#genome-file-as-input)
  - [**Protein and GFF file as input**](#protein-and-gff-file-as-input)
- [**Output**](#output)


# **About AnnotaPipeline**

Integrated tool to annotate hypothetical proteins developed by Laboratório de Bioinformática at Universidade Federal de Santa Catarina (Brazil).

AnnotaPipeline was tested in Unix-based systems. We strongly recommend to use Unix-based systems or Mac. AnnotaPipeline was not tested in Windows.

If you have questions, suggestions or difficulties regarding the pipeline, please do not hesitate to contact our team here on GitHub or by email: <labioinfo.genome@gmail.com>.


# **Analysis description**
> `AnnotaPipeline.py` checks all required parameters in `AnnotaPipeline.yaml` before execution

## **1. Gene prediction**

Run AUGUSTUS for gene prediction, using these required arguments:
  * strand
  * genemodel
  * species
  * protein
  * introns
  * start
  * stop
  * cds

You can also use these optional arguments:
  * hintsfile
  * extrinsicCfgFile
  * UTR

After gene prediction, sequences are "cleaned" based on minimal sequence size from `seq-cleaner` on `AnnotaPipeline.yaml`.

`.aa`, `.cdsexon` and `.codingseq` sequences are extracted from GFF output using `getAnnoFasta.pl` (AUGUSTUS script).

`.aa` sequences are used for subsequent analysis. `.codingseq` sequences are used for transcriptomics analysis (optional).


## **2. Similarity analysis**

Similarity analysis run through `blastp_parser.py`, which executes `blastp` with format 6 output (parsed after run) on predicted proteins (cleaned) against SwissProt:
  * qseqid
  * sseqid
  * sacc
  * bitscore
  * evalue
  * ppos
  * pident
  * qcovs
  * stitle
  > `evalue` and `max_target_seqs` given by user

This output is parsed to find annotations in the secondary database. The keyword list in `AnnotaPipeline.yaml` is used to exclude potential hypothetical annotations.

Hits are classified &ndash; as a potential annotation &ndash; if: (i) it doesn't have any words in the keyword list, and (ii) passed value thresholds for identity, positivity and coverage. All potential annotations &ndash; hits that passed all criteria &ndash; are present in `BASENAME_SwissProt_annotations.txt` for manual check.

Potential annotations are ranked by `bitscore` and the best hit is assigned to the corresponding protein. Proteins that aren't annotated based on SwissProt annotations, are separated in `BASENAME_BLASTp_AA_SwissProted.fasta` and `blastp` is rerun against the secondary database (selected in `AnnotaPipeline.yaml`).

Following the same process, `BASENAME_SpecifiedDB_annotations.txt` classifies potential annotations.

Annotation files:
  * `BASENAME_annotated_products.txt` contains all annotations
  * `BASENAME_hypothetical_products.txt` contains hypothetical proteins with no hits that passed all criteria(which doesn't have a single hit that passes all criteria) are present in `BASENAME_hypothetical_products.txt`
  * `BASENAME_no_hit_products.txt` contains proteins with no hits &ndash; against SwissProt and the secondary database &ndash; that are treated as hypothetical for subsequent analysis


## **3. Functional annotation**

**Functional analysis runs in two different ways:**
  * For annotated proteins:
    * InterProScan is used (with configured databases) to get ontology and IPR terms
  * For hypothetical proteins (which includes no_hit_products):
    * InterProScan, RPS-BLAST and HMMER are used to find hits of possible functions of predicted proteins

**Software arguments:**
  > Optional arguments can be given in `AnnotaPipeline.yaml` for InterProScan, HMMER and RPS-BLAST (not tested)
  * InterProScan &ndash; for annotated and hypothetical proteins &ndash; uses `-goterms` and `-iprlookup` arguments.
  * `hmmscan` runs with `--noali` argument and user values for `evalue` and `domE`. It also uses Pfam database.
  * RPS-BLAST runs with `evalue` and `max_target_seqs` arguments given in `AnnotaPipeline.yaml` and `-outfmt 6 "qseqid sseqid sacc bitscore evalue ppos pident qcovs stitle"`. RPS-BLAST uses CDD database.

**Parsing:**
  * `functional_annotation_parser.py` join outputs for both InterProScan runs in a single file called `InterProScan_Out_BASENAME.txt`. This file summarizes outputs for each predicted protein.
    * Coils, Gene3D and MobiDBLite databases are structural databases and are excluded from this output.
  * RPS-BLAST information for hypothetical proteins are summarized in `BASENAME_Grouped_Hypothetical_Information.txt` as it gives long descriptions. It may be helpful to find functional hints for proteins.
  * `info_parser.py` uses InterProScan results with parsed files from BLAST results to generate `All_annotated_products.txt`. This file joins gene annotation (from BLAST) with functional annotation (from InterProScan) using GOs and IPR values. This file is used to annotate sequences (aminoacid and nucleotide) in FASTA and GFF files.

## **4. Transcript quantification (Optional)**

Transcript quantification uses Kallisto and transcripts (`.codingseq` provided by AUGUSTUS) with RNA-seq data given by user.

Analysis starts with `kallisto index` and is followed by `kallisto quant`:
  > Both methods uses `bootstrap` value provided by user

  * 1 &ndash; paired-end data &ndash; using estimated average fragment length (`l`) and estimated standard deviation of fragment length (`s`)
    * Optional arguments for paired-end data if given by user in `AnnotaPipeline.yaml`

    OR

  * 2 &ndash; single-end data &ndash; using **REQUIRED** arguments estimated average fragment length (`l`) and estimated standard deviation of fragment length (`s`) given by user in `AnnotaPipeline.yaml`

**Parsing:**
  * `kallisto_parser.py` removes hits &ndash; by TPM value &ndash; below the `threshold` parameter from `proteomics` section. Possible thresholds are:
    * Median
    * Mean
    * Float value (user input)
  * Parsed output file is called `BASENAME_Transcript_Quantification.tsv` and contains a simplified result of Kallisto output (`abundance.tsv`) with target_id and TPM value

## **5. Peptide identification (Optional)**
> **WARNING:** Faild runs for Comet MS/MS don't crash AnnotaPipeline, so check AnnotaPipeline_Log.log to assure all spectrometry files produced outputs

Proteomics analysis uses Comet MS/MS with `comet.params` config given by user. In this file, our script overwrite values for following parameters:
  * `decoy_search = 1`
  * `output_pepxmlfile = 0`
  * `output_percolatorfile = 1`
  * `decoy_prefix = DECOY_`

Comet MS/MS runs with:
  * *Modified* `comet.params`
  * Annotated protein file
  * Path containing mass spectrometry files (`comet-spectrometry` param) and extension (`comet-ext` param)
    * Extension can be *mzXML, mzML, Thermo raw, mgf, and ms2 variants (cms2, bms2, ms2)*
  * Optional arguments `first` and `last` are used, if given, and overwrite cutoff values defined in `comet.params`

Comet MS/MS outputs (files with `.pin` extension) and raw files are moved to current directory.

Percolator runs for each `.pin` file with default parameters.

Percolator outputs are parsed by `percolator_parser.py` using `percolator-qvalue` (`AnnotaPipeline.yaml`). Raw files are maintained. `quantitative_proteomics` function uses parsed Percolator files to create `BASENAME_Total_Proteomics_Quantification.tsv`.

This output quantifies (all spectometry outputs):
  * Unique Peptide &ndash; number of unique peptides found across the entire dataset
  * Total Peptide &ndash; number of total peptides found across the entire dataset
  * Unique Spectrum &ndash; number of unique spectrum found across the entire dataset
  * Total Spectrum &ndash; number of unique spectrum found across the entire dataset

## **Last things last**

`summary_parser.py` integrates all outputs (in a single file summarizing all annotations found for each protein) from:
  * Prediction and similarity analysis
  * Functional annotation
  * Transcriptomic
  * Peptide identification (if used)


# **How to setup**

AnnotaPipeline requires the following software to run properly:
  * BLAST+ and RPS-BLAST (available at <https://ftp.ncbi.nih.gov/blast/executables/blast+/LATEST>)
  * InterProScan (available at <https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html>)
  * HMMER (available at <http://hmmer.org/download.html>)

You will also need to install [AUGUSTUS](https://doi.org/10.1093/bioinformatics/btr010) (available at <https://github.com/Gaius-Augustus/Augustus>) if you want to run this pipeline starting with gene/protein prediction.

Before executing, please modify the necessary fields in the configuration file (`AnnotaPipeline.yaml`).

## **Download databases**

**Required:**
  * SwissProt (available at <https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz>)
  * Pfam (available at <http://pfam.xfam.org>)
  * CDD (available at <https://www.ncbi.nlm.nih.gov/cdd>)

**Choose one secondary database:**
  * TrEMBL (available at <https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz>)
  * EuPathDB (available at <https://veupathdb.org/veupathdb/app>):
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

2. Create environment
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

2. Run `setup.py` (Scripts will be available at `$PATH`)

3. Install required softwares:
  * BLAST+ and RPS-BLAST (available at <https://ftp.ncbi.nih.gov/blast/executables/blast+/LATEST>)
  * InterProScan (available at <https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html>)
  * HMMER (available at <http://hmmer.org/download.html>)

4. Optional softwares
  * Kallisto (available at <https://pachterlab.github.io/kallisto/download.html>)
  * Comet MS/MS (available at <https://github.com/UWPR/Comet/releases/latest>)
    * Requires Percolator (available at <https://github.com/percolator/percolator>)

5. Download databases

6. Configure `AnnotaPipeline.yaml`

## **Setup `AnnotaPipeline.yaml`**
> **TIP:** If you already have InterProScan locally installed and configured, use it instead (as an alternative to conda installation &ndash; interpro vanilla)

We recommend using our example configuration file as a guide (`config_example.yaml`).


# **How to run**

`AnnotaPipeline.py` can run with three different options:
  * Protein file as input
  * Genome file as input
  * Protein and GFF file as input

## **Protein file as input**
```bash
AnnotaPipeline.py -p protein_sequences.fasta
```

This is the simplest execution of AnnotaPipeline.

The annotation process will begin with the submitted `protein_sequences.fasta` and a simplified version of `BASENAME_Annotated_GFF.gff` file at the end.

## **Genome file as input**
```bash
AnnotaPipeline.py -s genomic_data.fasta
```

This is the complete execution of AnnotaPipeline.

It will execute gene/protein prediction based on `genomic_data.fasta` utilizing AUGUSTUS and the predicted proteins will initiate the annotation process.

Given the prediction process, it is important to use a trained AUGUSTUS model for your species before executing AnnotaPipeline.

## **Protein and GFF file as input**
```bash
AnnotaPipeline.py -p protein_sequences.fasta -gff gff_file.gff
```

You can execute AnnotaPipeline &ndash; with this command line &ndash; if you already have a protein sequence file and a GFF3 file from previous AUGUSTUS predictions. The submitted `gff_file.gff` needs to be in GFF3 format.

The annotation process is the same as the first option, the difference being you will also have an annotated GFF3 output.


# **Output**
Depending on which option you decided to execute AnnotaPipeline, it will output three files (along with many others in their respective folders):
  * `All_Annotated_Products.txt` contains all unique sequence identifiers and their respective annotations (with functional annotations &ndash; when present).
  * ` Annota_BASENAME.fasta` contains all sequences and their annotations (with functional annotations &ndash; when present) in FASTA format.
  * `BASENAME_Annotated_GFF.gff` contains all sequences and their annotations (with functional annotations &ndash; when present) in GFF3 format.
  * `AnnotaPipeline_BASENAME_transcripts.fasta` contains nucleotide sequences for predicted proteins, with the same features present in the protein file.
  * `AnnotaPipeline_BASENAME_Summary.tsv` summarizes hits for each protein in similarity, functional, transcriptomics (if used) and proteomics analysis (if used).

Raw outputs are listed inside output folders:
  * `1_GenePrediction_BASENAME`           &ndash; AUGUSTUS files
  * `2_SimilarityAnalysis_BASENAME`       &ndash; BLAST+ analysis
  * `3_FunctionalAnnotation_BASENAME`     &ndash; InterProScan/HMMER/RPS-BLAST analysis
  * `4_TranscriptQuantification_BASENAME` &ndash; Kallisto analysis
  * `5_PeptideIdentification_BASENAME`    &ndash; Comet MS/MS and Percolator analysis

The output folders and files will be located in the same folder you executed the pipeline.