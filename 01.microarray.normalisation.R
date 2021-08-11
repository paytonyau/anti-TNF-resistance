BiocManager::install("affy")
library(affy)
library(tcltk)

## Download the corresponding raw data from GEO database
# Apply choose.dir  to select a directory
dir <- tk_choose.dir(caption = "Select folder")
# List CEL files and save the variables
cel.files <- list.files(path = dir, pattern = ".+\\.cel$", ignore.case = TRUE, full.names = TRUE, recursive = TRUE)
# View the file names
basename(cel.files)
data.raw <- ReadAffy(filenames = cel.files)

# background processing using rma function
# normalisation processing uses quantile method 
# summary method uses medianpolish
eset.rma <- rma(data.raw)

# format the matrix
norm.rma <- data.frame(eset.rma)
norm.rma <- t(norm.rma)
norm.rma <- data.frame(norm.rma)

# Back up the matrix
write.csv(norm.rma, file = "_rma.csv") 
# Probes with capital letter "X" on the first letter are required to remove 
# in order to let the downstram program reconise the probes without "X"
df2 <- read.csv("_rma.csv")

# Annotation 
BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
## Change it based on microarray model
library(hgu133plus2.db) ## For GSE16879/GSE23597/GSE52746
library(hugene10stv1probe) ## For GSE92415
## the annotation for GSE73661 first used the annotation matrix from GSE database
## and update the gene names using HGNChelper R package

## From Colname to the first column
K <- df2$X
K <- as.vector(K)

## PROBE to Gene
SYMBOL <- mapIds(hgu133plus2.db, keys=K, column= "SYMBOL", keytype="PROBEID", multiVals="first")
df2$X <- SYMBOL

## omit NA values
df3 <- na.omit(df2, cols = df2$X)

#if multiple probe sets map to a gene, select the one with median
C = aggregate(df3, by = list(df3$X), FUN = median)
rownames(C) <- C$Group.1
C$X <- NULL
C$Group.1 <- NULL
write.csv(C, file = "_rma_median.csv")
