library(edgeR)
clic <- read.csv("GSE145918_Clin.csv", header = T, row.names = 1)
counts <- read.delim("GSE145918_raw_counts_2.txt", row.names = 1)
d0 <- DGEList(counts, group = clic$group)
d0$samples

## Filter low-expressed genes
d <- filterByExpr(d0)
d0 <- d0[d, , keep.lib.sizes=FALSE]
dim(d0) # number of genes left

##
d <- d0

#### 2. Preprocessing
## calcNormFactors doesn't normalize the data, 
## it just calculates normalization factors for use downstream.
d0 <- calcNormFactors(d0) # TMM, RLE or upperquartile, default TMM
CPM <- cpm(d0, prior.count=3) # log2 transformation
logCPM <- log2(CPM + 1)
write.csv(logCPM, "log2_CPM_Matrix.csv")

mod = model.matrix(~as.factor(group), data=clic)
mod0 = model.matrix(~1,data=clic)

library(sva)
svobj = svaseq(as.matrix(d$counts), mod, mod0)
modSv = cbind(mod,svobj$sv)

## 3. Voom transformation and calculation of variance weights
mm <- model.matrix(~0+group, data = clic)
mm = as.data.frame(mm)
modSv = cbind(mm,svobj$sv)

### change 1, 2 to A- B.....
colnames(modSv)[7:8] = c("A", "B")

## Since we need to make comparisons both within and between subjects, 
## it is necessary to treat Patient as a random effect

## https://stats.stackexchange.com/questions/160255/voom-mean-variance-trend-plot-how-to-interpret-the-plot
y <- voom(d, modSv, plot = T)

## 4. Fitting linear models in limma
fit <- lmFit(y, modSv)
head(coef(fit))

contr <- makeContrasts( bn_ctr_tnf = groupbn_ctr - groupbn_tnf, 
                        bn_ctr_lps = groupbn_ctr - groupbn_lps,
                        bn_tnf_lps = groupbn_tnf - groupbn_lps,
                        dn_ctr_tnf = groupdn_ctr - groupdn_tnf,
                        dn_ctr_lps = groupdn_ctr - groupdn_lps,
                        dn_tnf_lps = groupdn_tnf - groupdn_lps,
                       levels = colnames(coef(fit)))
contr

# Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)
tmp2 <- eBayes(tmp)

top.table <- topTable(tmp2, coef="dn_tnf_lps", sort.by = "P", n = Inf)

length(which(top.table$adj.P.Val < 0.05)) # number of DE genes
write.table(top.table, file = "dn_tnf_lps.txt", sep = "\t", quote = F)


##############################################################
library(EnhancedVolcano)
## Adjusted P-Value Plot
Volcan1<- EnhancedVolcano(top.table,
                          lab = row.names(top.table),
                          x = 'logFC',
                          y = 'P.Value',
                          title = "Refractory vs Onset",
                          subtitle = "GvHD patients",
                          ylab = bquote(~-Log[10]~italic(P)~value),
                          xlim = c(-2, 2),
                          ylim = c(0, 5),
                          shape = 16,
                          pointSize = 3.5,
                          colAlpha = 4/5,
                          # col=c('black', 'green', 'gray', 'red3'),
                          legendPosition = "bottom",
                          legendLabSize = 14,
                          drawConnectors = FALSE,
                          pCutoff = 0.01,
                          FCcutoff = 0.75)

pdf("rplot6.pdf", width = 10, height = 10) 
Volcan1  # + ggplot2::scale_x_continuous(
# breaks=seq(-2, 2, 1))
dev.off()
