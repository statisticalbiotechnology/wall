if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install()
BiocManager::install(c("doParallel", "BiocParallel", "DESeq2", "devtools", "edgeR"))
install.packages('devtools', repos = "http://cran.us.r-project.org")

#devtools::install_git("https://git.bioconductor.org/packages/edgeR")
library("devtools")
install_github("statOmics/zingeR")

install_github('xzhoulab/iDEA', dependencies=TRUE)
