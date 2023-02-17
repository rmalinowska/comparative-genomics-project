# Project for comparative genomics course

Project aims to perform analysis of given species proteomes and create proposed phylogenetic trees using different methods like consensus tree and supertree. 

### How to run
Required packages: Python: Biopython, numpy, matplotlib, requests, copy, sys, re, json, os; R: TreeTools, ape.
Required programms: muscle, mmseqs2, duptree.

#### Data
There are two options of providing input data to the script:

1. file with species names separated with new line, for example:

Aspergillus nidulans

Candida parapsilosis

Magnaporthe grisea

Debaryomyces hansenii

Candida tenuis

2. file with species names and proteome IDs from InterPro datapase separated with new line, for example:

Aspergillus nidulans,UP000000560

Candida parapsilosis,UP000005221

Schizosaccharomyces japonicus,UP000001744

First option is good for the case, where we have a lot of species and it would be tedious to search for all proteome IDs, but in some cases it might result with inproper results, cosidering that a wrong proteome may be downloaded from database (it downloads the first result from browsing proteomes by species name). In this case it is recommended to validate results (for example check if size of downloaded proteome is not suspiciously small).

Second option makes sure that exact proteomes will be downloaded, but it requires providing them in the first place.

Using first option:

In the beginning of the script you have to set variable SPECIES to True and provide name of the file and leave variable IDS set to False.

Using second option:

In the beginning of the script you have to set variable IDS to True and provide name of the file and leave variable SPECIES set to False.

#### The rest of parameters:

In the beginning of the script you have to also provide multiple other parameters:

working_dir = *str: path to directory where you store input file and where results will appear*

merged_fasta_output_filename = *str: file name which will merge all proteomes*

path_to_mmseq = *str: path to mmseqs2*

IDENT = *float: minimal identity threshold for clustering*

COVER = *float: coverage percentage for clustering*

COV_MODE = *int: coverage mode for clustering*

PREF = *str: prefix for clustering output files*

path_to_muscle = *str: path to muscle*

TREE_CONSTRUCTION_METHOD = *str: method for buliding trees in Phylo.TreeConstruction, for example 'nj' or 'mp'*

DISTANCE_METHOD = *str: method for calculating distance between sequences for tree construction, for example: 'blosum62', 'pam70'*

PATH_RSCRIPT = *str: path to Rscript in orter to run R script; while using linux operating system, it is probably sufficient to set it to 'Rscript'*

PATH_DUPTREE = *path to duptree*

BOOTSTRAP_REP = *int: number of bootstrap samples*

To run the script we have to store following files in the working directory: projekt_gp.py, consensus.r, download_proteome.py and run projekt_gp.py file.
