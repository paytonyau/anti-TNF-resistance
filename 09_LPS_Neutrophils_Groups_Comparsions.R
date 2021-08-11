#### https://www.datanovia.com/en/blog/how-to-perform-t-test-for-multiple-variables-in-r-pairwise-group-comparisons/
### Load required R packages
library(tidyverse)
library(rstatix)
library(ggpubr)

#### Pre-calculated the mean expression level (across samples) of genes matching indicated GO term (Figure 3D-3E) from the IBD patients
mydata.long <- read.csv("Chemoteaix2.csv")

## stat.test <- mydata.long %>%
##  group_by(Group) %>%
##   tukey_hsd(Avg.Value ~ Treatment)

stat.test <- mydata.long %>%
  group_by(Group) %>%
  wilcox_test(Avg.Value ~ Treatment) %>%
  add_significance("p")
stat.test

## Remove unnecessary columns and display the outputs
stat.test %>% select(-.y., -statistic, -df)
stat.test

## Visualisation
pdf("LPS-Neutrophils.pdf",7 ,6)
myplot <- ggbarplot(
  mydata.long, x = "Treatment", y = "Avg.Value",
  facet.by = "Group", add = "mean_se",size = 0.5, 
  linetype = 1,
  # color= "Treatment", 
  ylab = "Mean value of gene expression (Log2 + 1)", 
  xlab = FALSE,
  fill = "Treatment",
  ggtheme = theme_pubr(border = TRUE)
) + 
  facet_wrap(~Group) + #  theme_minimal() + # + scale_fill_brewer(palette="Dark2")
  theme_classic() + 
  theme(text = element_text(size=15, colour = "black"), 
                            axis.ticks = element_line(colour = "black", size = 1),
                            axis.line = element_line(colour = 'black', size = 1),
                            axis.text.x = element_text(angle=0, hjust=0.5, colour = "black",
                                                       size = 13, face="bold"),
                            axis.text.y = element_text(angle=0, hjust=0.5, colour = "black",
                                                       size = 13, face="bold"),
                            axis.title.y = element_text(color="black", size=15,face="bold"),
                            legend.position = "none") + 
  scale_y_continuous(limits=c(0, 11), breaks = c(0, 2.5, 5, 7.5, 10))

## Add statistical test p-values
stat.test <- stat.test %>% 
  add_xy_position(x = "Treatment")

myplot + stat_pvalue_manual(stat.test, label = "p.signif",
                            y.position = stat.test$y.position + c(0.2, 0.7, 1.3))
dev.off()
