# COMS4761_Project
Final project for COMS4761 - Computational Genomics

In order to run this project, a python code is first used to generate the
initial results by performing the alignments and then manipulating these
data into formats that will be used for analysis by R scripts. The alignment
files are all already given here, so those steps can be skipped when testing
this project. The R scripts are broken up by the type of analysis being run
and generates output files appropriately.

This project does not require a powerful machine to run and can be executed
on most machines as data manipulation has been converted to summarized 
data in .txt formats. 


The Files Available on the Git:

>>> ISH_pTRAP
Contains a folder that contains exploratory scatter plots and exp_plot
(expression plots). These are such simple scatter plots based on in-situ 
hybridization and PhosphoTrap metrics available in an attempt to find a 
linear relationship between the two data. 

plot_script.R makes differential expression plots from the sleuth object
generated in the above code.

IDR.R is used to generate the t vs. psi plots of in-situ hybridization data
against PhosphoTrap data

IDR_PTrap is to generate the IDR results of PhosphoTrap replicates. T

>>> alignment_results
This folder contains the alignment results of all the kallisto alignments performed
on the RNAseq data that we had available. As such, the folder contains sub folders
that are the names of the samples (signifying the type of data and stimulus) which
contains the output of the alignment itself. It also CR_metadata.txt, which is just 
a text file of all the folders associated with control data. This was just used to 
automate some of the data manipulation that had to be done with these files. T

>>> conversation_tables
This folder contains the txt files of the conversions available between various 
methods of denoting genes/transcript. The naming convension of the file signifies 
the conversion stored in the text file. Conversion_url.txt is simply the link to 
the database by which all the conversion files were taken.

>>> ranked_metrics
This folder contains csv files where the values of interest are listed in the order
of ranks with their associated gene/transcript ID. This was used for the purposes 
of IDR analysis.

>>> sleuth_results
As the name suggests, this is a folder containing results of the sleuth pipeline and
subsequent analysis. As such, it contains the volcano plots and DE rankings of the genes
based on beta value generated from a Wald Test.

>>> Python_Code:

This is the Jupyter Notebook where a lot of the pre-processing of data files
occurred. Much of the statistical analysis was done using R. As such, this
notebook contains functions for manipulating files and data. Each function 
themselves has been commented with what it does along with the input para-
meters and return values. To summarize, the general workflow is:

1) index_transcriptome: generates a indexed reference file for kallisto 
alignment. The FASTA reference used and the index file are too large to be
included, but the files have been taken from:

https://useast.ensembl.org/info/data/ftp/index.html

The FASTA file used is the cDNA of mouse, specifically:

Mus_musculus.GRCm38.cdna.all.fa.gz

2) mass_alignments: generates the kallisto alignments of all the RNAseq data
against the indexed reference generated above. This was done to automate the
process and was done with a bootstrap value of 100. Running all the alignments
takes approximately 1 hour for each sample (30+ hours total). However, 
the aligned files are already included for convenience under the alignment_result
folder that is included. 

3) conversion_dict: makes python dictionaries of conversion. These conversion files
are all located in the conversion_table folder.

4) Convert all the data into formats that can be manipulated in Python

phos_meta: returns a dictionary of averaged data from all the control samples using
the CR_metadata file available.

ish_data: returns a dictionary of the in-situ hybridization data from the Allen B
Brain Atlas 

ish_phos_data: compiles all the data from the phos_meta dict and the ish_data dict
using the conversion_dict to map between the two different types of data

convert_to_csv: takes the data outputted from ish_phos_data and puts it into csv 
format for manipulation in R since Python code and R code are not linked with 
another script. 

5) The remaining part of the notebook is doing analysis of the data and constructing
a summary of all the data, including ranks for further analysis

>>> ish_phos_merge.csv
csv file generated from the Python pipeline that contains all of the data from
the Allen Brain Atlas and PhosphoTrap data mapped to each other


>>> sleuth_script.R
sleuth_script.R is an implementation of the sleuth pipeline that generates
sleuth objects


