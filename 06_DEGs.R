# install limma package
install::packages("limma")
library("limma")

# import the matrix
Z2 = read.csv("matrix.csv", header = TRUE, rowname = 1)

# create design matrix for the model
design_palmieri <- model.matrix(~ 0 + Z2$A) # Group

# set the column names of the design matrix
colnames(design_palmieri)[1:2] <- c("IFN_Pre_NR", "IFN_Pre_R")

# set the row names of the design matrix
rownames(design_palmieri) <- Z2$GSM

# create contrast matrix
contrast_matrix <- makeContrasts(IFN_Pre_NR-IFN_Pre_R, levels = Z2$A)

# fit the linear model
palmieri_fit <- eBayes(contrasts.fit(lmFit(combat_edata2,design = design_palmieri),contrast_matrix))

# create a table of the top results
table <- topTable(palmieri_fit, number = Inf)

# install and load EnhancedVolcano package
install::packages("EnhancedVolcano")
library(EnhancedVolcano)

## create an adjusted p-value volcano plot
Volcan1<- EnhancedVolcano(table,
                          lab = row.names(table),
                          x = 'logFC',
                          y = 'adj.P.Val',
                          title = "Responder vs Non-Responder",
                          subtitle = "Pre TNF-Alpha Treatment",
                          ylab = bquote(~-Log[10]~adjusted~italic(P)~value),
                          #xlim = c(-2, 2),
                          ylim = c(0, 10),
                          shape = 16,
                          pointSize = 3.5,
                          colAlpha = 4/5,
                          # col=c('black', 'green', 'gray', 'red3'),
                          legendPosition = "bottom",
                          legendLabSize = 14,
                          drawConnectors = TRUE,
                          pCutoff = 0.05,
                          FCcutoff = 0.75)

# save the volcano plot as a pdf file
pdf("Volcano.pdf", width = 10, height = 10) 
Volcan1  # + ggplot2::scale_x_continuous(
  # breaks=seq(-2, 2, 1))
dev.off()

# create a subset of the data where logFC is greater than or equal to 0.75 and adj.P.Val is less than or equal to 0.05,
# or where logFC is less than or equal to -0.75 and adj.P.Val is less than or equal to 0.05
sub_group <- subset(A, logFC >= 0.75 & adj.P.Val <= 0.05 | 
                      logFC <= -0.75 & adj.P.Val <= 0.05)

# add the Gene column to both sub_group and combat_edata
sub_group$Gene <- row.names(sub_group)
combat_edata$Gene <- row.names(combat_edata)

# merge the sub_group and combat_edata data frames on the Gene column
exprset <- merge(sub_group, combat_edata, by.x="Gene", by.y="Gene")

# set the row names of the exprset data frame
rownames(exprset) <- exprset$Gene

exprset <- exprset[,8:ncol(exprset)]

#############
## Creating Heatmap using pheatmap package
library(pheatmap)
library(RColorBrewer)

# Create a z_score loop function
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# z-score calculation
data_subset <- t(apply(exprset, 1, cal_z_score))

# some data re-formatting
data_subset <- as.data.frame(data_subset)
data_subset <- t(data_subset)
data_subset <- as.data.frame(data_subset)

annotation_col = as.data.frame(Z2$A)
colnames(annotation_col) <- c("Status")

library(dplyr)

annotation_col<- annotation_col %>% 
  mutate(Status = ifelse(as.character(Status) == "IFN_Pre_R", 
                          "Resistance", as.character(Status)))

annotation_col<- annotation_col %>% 
  mutate(Status = ifelse(as.character(Status) == "IFN_Pre_NR", 
                         "Non-Resistance", as.character(Status)))
rownames(annotation_col) <- Z2$GSM

# Plot the heatmap using pdf format
pdf("Heatmap.pdf",width = 25, height = 12.5)
pheatmap(data_subset, color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         clustering_method = "average", annotation_col = annotation_col)
dev.off()

###################
## PCA analysis
library(FactoMineR)
library(factoextra)

## Prepare the data
data = as.data.frame(t(data_subset))
data <- cbind(data, Z = as.factor(Z2$A))

pca <- PCA(data[,-ncol(data)], graph = FALSE)
eig.val <- get_eigenvalue(pca)
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 100))

# plot PCA
pdf("pca_before.pdf",width = 8, height = 6)
fviz_pca_ind(pca,
             geom.ind = "point", 
             col.ind = annotation_col$Status, 
             palette = c("dodgerblue4", "burlywood4", "#56B4E9"), 
             addEllipses = TRUE, 
             legend.title = "TNF-Alpha Treatment",
             title = '')
dev.off()
