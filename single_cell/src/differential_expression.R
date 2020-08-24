#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install(c("doParallel", "BiocParallel", "RCurl", "DESeq2", "zingeR", "devtools"))

library("DESeq2")
#library("DFP")
library("doParallel")
library("BiocParallel")
#BiocManager::install("BiocParallel")
#library("Biobase")
#BiocManager::install("DESeq2")

#library("zinbwave")
#library("scales")
#library("limma")
#library("MAST")
#library("edgeR")
library("devtools")


#devtools::install_git("https://git.bioconductor.org/packages/edgeR")
#install_github("statOmics/zingeR")
library("zingeR")

#BiocManager::install("zingeR")
#setwd("C:/Users/patri/Documents/GitHub/wall/single_cell/data/finished_data_files")
#NCORES <- 1
#registerDoParallel(NCORES)
#register(DoparParam())


full_df <- read.csv("full_df.csv", sep = ",", header=TRUE, row.names=1)
matrix_df <- data.matrix(full_df)



pheno_data <- read.csv("../pheno_df.csv", sep=",", row.names = 1, header = TRUE)
metaData <- data.frame(labelDescription=c("PatientType"))
annData <- AnnotatedDataFrame(data=pheno_data, varMetadata = metaData)
pheno_data <- as(pheno_data, "AnnotatedDataFrame")


eset <- ExpressionSet(matrix_df, pheno_data)
counts = exprs(eset)
cellType =pData(eset)[,"phenotype"]

colData <- data.frame(cellType = cellType)
design <- model.matrix(~ cellType)
dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = design)
weights <- zingeR::zeroWeightsLS(counts = counts, design = design, maxit = 500, normalization = "DESeq2_poscounts", colData = colData, designFormula = ~cellType, verbose = FALSE)
assays(dse)[["weights"]] <- weights
dse <- DESeq2::estimateSizeFactors(dse, type="poscounts")
dse <- estimateDispersions(dse)
dse <- nbinomWaldTest(dse, modelMatrixType="standard", betaPrior=TRUE, useT=TRUE, df=rowSums(weights)-13, minmu=1e-3)


install.packages("NCmisc")
library("NCmisc")
list.functions.in.file("data_perparation.R", alphabetic = TRUE)






