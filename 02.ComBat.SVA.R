# Load the BiocManager package if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install the sva package
BiocManager::install("sva")

# Load the sva package
library(sva)

# Input the combined clinical information
sample.list <- read.csv("Sample_list_filtered.csv", header = TRUE)
batch <- sample.list$Batch

# Input the merged matrix from five different datasets
df <- read.csv("merged_matrix.csv")

# Apply the null model
modcombat <- model.matrix(~1,data=df)

# Adjust for batch with combat function
combat_edata <- ComBat(dat = as.matrix(df), batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)

# Transpose the data and convert to a data frame
combat_edata <- as.data.frame(t(combat_edata))

# Save the batch effect corrected matrix
write.csv(combat_edata, "merged_batch_effect_corrected_matrix.csv")

# Plot figures
pdf("Batch_effect_After_correction.pdf", width = 30, height = 10)
boxplot(combat_edata)
dev.off()
pdf("Batch_effect_Before_correction.pdf", width = 30, height = 10)
boxplot(df)
dev.off()
