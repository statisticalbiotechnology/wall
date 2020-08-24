setwd("/mnt/c/Users/patri/Documents/GitHub/wall/single_cell/data/finished_data_files")
library(iDEA)



summary_stats <- read.csv("lfc_df.csv", sep = ",", header=TRUE, row.names=1)
annotation_data <- read.csv("annotation_reactome.csv", sep= ",", header = TRUE, row.names = 1)
data(summary_stats)
head(summary_stats)



idea <- CreateiDEAObject(summary_stats, annotation_data, max_var_beta = 100, min_precent_annot = 0.00025, num_core=10)

head(idea@summary)
head(idea@annotation[[1]])

idea <- iDEA.fit(idea,
                 fit_noGS=FALSE,
                 init_beta=NULL, 
                 init_tau=c(-2,0.5),
                 min_degene=1,
                 em_iter=15,
                 mcmc_iter=1000, 
                 fit.tol=1e-5,
                 modelVariant = F,
                 verbose=TRUE)

idea@de
head(idea@de)

idea <- iDEA.louis(idea)

idea@gsea
head(idea@gsea)

write.csv(idea@gsea,"idea_pathways.csv", row.names = FALSE)
idea <- iDEA.BMA(idea)
head(idea@BMA_pip)


idea_variant <- iDEA.fit(idea,modelVariant = T)


idea.null <- CreateiDEAObject(summary_stats, annotation_data[,c(1:10)], num_core=10)
idea.null <- iDEA.fit.null(idea.null,numPermute = 10) ## 
idea.null <- iDEA.louis(idea.null) 

head(idea.null@gsea)

df.FDR <- iDEA.FDR(idea,idea.null,numPermute = 10)
head(df.FDR)
write.csv(df.FDR, "idea_pathway_fdr.csv", row.names=FALSE)
