if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("sva")

library(sva)

sample.list = read.csv("Sample_list_filtered.csv", header = T)
batch = sample.list$Batch

df <- read.csv("....csv")

### null model
modcombat = model.matrix(~1,data=df)
### Adjust for batch with combat function ## https://www.biostars.org/p/142914/
combat_edata = ComBat(dat=as.matrix(df), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
combat_edata <- as.data.frame(t(combat_edata))
write.csv(combat_edata, "A.csv")

## Plot figures
pdf("Batch_effect_After_correction.pdf", 30, 10)
boxplot(combat_edata)
dev.off()
pdf("Batch_effect_Before_correction.pdf", 30, 10)
boxplot(df)
dev.off()