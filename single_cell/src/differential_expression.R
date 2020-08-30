library("DESeq2")
library("doParallel")
library("BiocParallel")
library("zingeR")


full_df <- read.csv("../data/full_df.csv", sep = ",", header=TRUE, row.names=1)
matrix_df <- data.matrix(full_df)


pheno_data <- read.csv("../data/pheno_df.csv", sep=",", row.names = 1, header = TRUE)
metaData <- data.frame(labelDescription=c("phenotype"))
annData <- AnnotatedDataFrame(data=pheno_data, varMetadata = metaData)
pheno_data <- as(pheno_data, "AnnotatedDataFrame")


eset <- ExpressionSet(matrix_df, pheno_data)
counts = exprs(eset)
cellType = pData(eset)[,"phenotype"]

colData <- data.frame(cellType = cellType)
design <- model.matrix(~cellType)
dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~cellType)
weights <- zingeR::zeroWeightsLS(counts = counts, design = design, maxit = 500, normalization = "DESeq2_poscounts", colData = colData, designFormula = ~cellType, verbose = TRUE)
assays(dse)[["weights"]] <- weights
dse <- DESeq2::estimateSizeFactors(dse, type="poscounts")
dse <- estimateDispersions(dse)
dse <- nbinomWaldTest(dse, modelMatrixType="standard", betaPrior=TRUE, useT=TRUE, df=rowSums(weights)-13, minmu=1e-3)
resultsNames(dse)

res <- results(dse)
res <- res[c("log2FoldChange", "lfcSE")]

write.csv(as.data.frame(res), file="../data/lfc_df.csv")
