# The data and code for paper "Widespread transposable element dysregulation in human aging brains with Alzheimer’s disease"


## Data Description
This document contains scripts for analyzing transposable element (TE) dysregulation of quantitative trait loci (QTLs) in human aging brains with Alzheimer’s disease. It contains four scripts


### 1. Transposable element expression quantification using bulk RNA-seq data

#### Software requirement:
* SQuIRE (v0.9.9.9a-beta): https://github.com/wyang17/SQuIRE
* Telescope (v1.0.3): https://github.com/mlbendall/telescope
* TEtranscripts (v2.2.3): https://github.com/mhammell-laboratory/TEtranscripts

#### Files
* `SQuIRE` folder includes codes used for quantifying transposable element expression by SQuIRE
* `Telescope` folder includes codes used for quantifying transposable element expression by Telescope
* `TEtranscripts` folder includes codes used for quantifying transposable element expression by TEtranscripts


### 2. TE expression quantitative trait loci (teQTL)

#### Software requirement:
* MatrixEQTL: https://github.com/andreyshabalin/MatrixEQTL

#### Files
* `Mayo_teQTL.txt.gz`: teQTLs identified from the 152 temporal cortex (TCX) brains in Mayo brain biobank
* `README.txt`: detailed description of Mayo_teQTL.txt.gz file
* `teQTL_pepline_using_Matrix_eQTL.R` is the code used for teQTL identification.


### 3. Colocalization analysis

#### Software requirement:
* coloc (5.1.0.1): https://github.com/chr1swallace/coloc


* Codes used for transposable element expression at single-cell level have been uploaded to below link: https://github.com/ChengF-Lab/alzTE-HHV/tree/main/02_scTE_analysis
