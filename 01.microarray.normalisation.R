BiocManager::install("affy")
library(affy)
library(tcltk)

## Before starting the process
## Download the corresponding raw data from GEO database

# Select a directory using tk_choose.dir
dir <- tk_choose.dir(caption = "Select folder")

# Get the list of CEL files in the selected directory
cel.files <- list.files(path = dir, pattern = ".+\\.cel$", ignore.case = TRUE, full.names = TRUE, recursive = TRUE)

# Identify the file names
basename(cel.files)

# Read the data into a raw data object
data.raw <- ReadAffy(filenames = cel.files)

# Normalize and summarise the raw data using rma method with quantile normalization
eset.rma <- rma(data.raw)

# Transpose and format the resulting matrix as a data frame
norm.rma <- data.frame(eset.rma)
norm.rma <- t(norm.rma)

# Write the matrix to a CSV file
norm.rma <- data.frame(norm.rma)

# Back up the matrix
write.csv(norm.rma, file = "_rma.csv") 
# Probes with capital letter "X" on the first letter are required to remove 
# in order to let the downstram program reconise the probes without "X"

# Read the CSV file into a data frame
df2 <- read.csv("_rma.csv")

# Annotation: convert the probe IDs to gene symbols using mapIds and the hgu133plus2.db or hugene10stv1probe library
BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
## Change it based on the microarray model

## (1) For GSE16879/GSE23597/GSE52746
library(hgu133plus2.db)
## (2) For GSE92415
library(hugene10stv1probe)

## (3) For GSE73661 
# the annotation first used the annotation matrix from the GSE database
## and update the gene names using HGNChelper R package

## From Colname to the first column
K <- df2$X
K <- as.vector(K)

## PROBE to Gene
SYMBOL <- mapIds(hgu133plus2.db, keys=K, column= "SYMBOL", keytype="PROBEID", multiVals="first")
df2$X <- SYMBOL

# Omit rows with NA values
df3 <- na.omit(df2, cols = df2$X)

# Aggregate the data frame by gene, selecting the median value for each gene
C = aggregate(df3, by = list(df3$X), FUN = median)
rownames(C) <- C$Group.1
C$X <- NULL
C$Group.1 <- NULL
write.csv(C, file = "_rma_median.csv")
