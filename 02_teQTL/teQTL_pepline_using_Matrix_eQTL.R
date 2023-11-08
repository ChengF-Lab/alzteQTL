#!/usr/bin/Rscript

args <- commandArgs(TRUE)
geno_dir <- args[1]
rpkm_dir <- args[2]
eqtl_dir <- args[3]
chr <- args[4]

library(MatrixEQTL)

## Settings
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR

# Genotype file name
SNP_file_name = paste0(geno_dir,"/",chr,"/Mayo_",chr,"_geno.info")

# Gene expression file name
expression_file_name = paste0(rpkm_dir,"/",chr,"/Mayo_",chr,"_transcript.txt")

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste0(eqtl_dir,"/",chr,"/Mayo_",chr,"_covaria.txt")

# Output file name
output_file_name = paste0(eqtl_dir,"/",chr,"/Mayo_",chr,"_trans_eQTL_1Mb.txt")
output_file_name.cis = paste0(eqtl_dir,"/",chr,"/Mayo_",chr,"_cis_eQTL_1Mb.txt")

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric()

## Load genotype data
snps = SlicedData$new()
snps$fileDelimiter = "\t"              # the TAB character
snps$fileOmitCharacters = "NA"         # denote missing values;
snps$fileSkipRows = 1                  # one row of column labels
snps$fileSkipColumns = 1               # one column of row labels
snps$fileSliceSize = 2000              # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)

## Load gene expression data
gene = SlicedData$new() 
gene$fileDelimiter = "\t"              # the TAB character
gene$fileOmitCharacters = "NA"         # denote missing values;
gene$fileSkipRows = 1                  # one row of column labels
gene$fileSkipColumns = 1               # one column of row labels
gene$fileSliceSize = 2000              # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)

## Load covariates
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"              # the TAB character
cvrt$fileOmitCharacters = "NA"         # denote missing values;
cvrt$fileSkipRows = 1                  # one row of column labels
cvrt$fileSkipColumns = 1               # one column of row labels
if(length(covariates_file_name)>0){
	    cvrt$LoadFile(covariates_file_name);
}

## Load snps position files (data.frame format)
snpspos <- read.table(paste0(geno_dir,"/",chr,"/Mayo_",chr,"_snp.pos"),header = T)
snpspos <- as.data.frame(snpspos)

## Load expresison position files (data.frame format)
genepos <- read.table(paste0(rpkm_dir,"/",chr,"/Mayo_",chr,"_transcript.pos"),header = T)
genepos <- as.data.frame(genepos)

## Run the analysis

me = Matrix_eQTL_main(
      snps = snps, 
	  gene = gene, 
	  cvrt = cvrt, 
	  output_file_name = output_file_name, 
	  pvOutputThreshold = 0,
	  useModel = modelLINEAR, 
	  errorCovariance = numeric(), 
	  verbose = TRUE, 
	  output_file_name.cis = output_file_name.cis, 
	  pvOutputThreshold.cis = 1e-1,
	  snpspos = snpspos, 
	  genepos = genepos,
	  cisDist = 1e6,
	  pvalue.hist = "qqplot",
	  min.pv.by.genesnp = FALSE,
	  noFDRsaveMemory = FALSE)


#cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
#cat('Detected eQTLs:', '\n');


