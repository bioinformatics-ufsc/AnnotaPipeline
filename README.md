# AnnotaPipeline
Integrated tool to annotate hypothetical proteins developed by the Laboratorio de Bioinformatica at the Universidade Federal de Santa Catarina, in Brazil.  
If you have questions about the pipeline or difficulties with its execution, please do not hesitate to contact our team here on GitHub or to contact me by email: Guilherme Maia <labioinfo.genome@gmail.com>

# How to setup
The AnnotaPipeline requires the following software to run properly:  
BLAST and RPSblast (available at https://ftp.ncbi.nih.gov/blast/executables/blast+/LATEST/)  
InterProScan (available at https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html)  
HMMER (available at http://hmmer.org/download.html)  
You will also need to install Augustus (available at https://github.com/Gaius-Augustus/Augustus) if you want to run this pipeline starting from gene and protein prediction.   
Before executing, please, modify the necessary fields in the AnnotaPipeline.config file.  

# How to run
To initiate the pipeline, you can run the ```AnnotaPipeline.py``` script with three different options:  

```python3 AnnotaPipeline.py -p protein_sequences.fasta```  
This is the most simple version of AnnotaPipeline execution.  
The annotation process will begin with the submitted ```protein_sequences.fasta``` and the only missing output will be the ```(basename)_Annotated_GFF.gff``` file at the end.  

```python3 AnnotaPipeline.py -s genomic_data.fasta```  
This is the most complete version of AnnotaPipeline execution.  
It will execute gene and protein prediction on the ```genomic_data.fasta``` utilizing Augustus, and then utilize the predicted proteins to initiate the annotation process.  
Because of this first step, it is important that you have trained Augustus to your particular species before executing AnnotaPipeline.  

# Output
Depending on which option you decided to execute the AnnotaPipeline, it will output three files (among many others in their respective folders):  
```All_Annotated_Products.txt``` which contains all unique sequence identifiers and their respective annotations (with functional annotations, when present).  
``` Annota_(basename).fasta``` which contains all the sequences and their annotations (with functional annotations, when present) in FASTA format.  
```(basename)_Annotated_GFF.gff``` which contains all the sequences and their annotations (com anotações funcionais, quando presente) in GFF3 format.  
The output folders and files will be located in the same folder you executed the pipeline in.  
