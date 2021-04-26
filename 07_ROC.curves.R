library(tidyverse)
library(pROC)

ROCStatFunc <- function(dat, group, var,retype = c("threshold", "specificity", "sensitivity"),
                        auc = T,youden = T, digit = 3){
  subgroup <- levels(as.factor(dat[[group]]))
  subgroup1 <- paste0(subgroup[2], " vs ", subgroup[1])
  rocmodel <- roc(dat[[group]], dat[[var]])
  other <- coords(rocmodel, "b", ret = retype)
  other <- round(other, digit)
  if(auc == T){
    auc <- round(ci.auc(rocmodel),digit)
    auc <- paste0(auc[2],"(",auc[1],"-",auc[3],")")
    if(youden == T){
      abc <- coords(rocmodel, "b", ret = c("specificity", "sensitivity"))
      youdenres <- abc[1] + abc[2] - 1
      youdenres <- round(youdenres, digit)
      result <- c(group, subgroup1, auc, other, youdenres)
      names(result) <- c("group", "subgroup","auc(95%CI)", retype, "youden")
    }else{
      result <- c(group, subgroup1, auc, other)
      names(result) <- c("group", "subgroup", "auc(95%CI)", retype)
    }
  }else{
    if(youden == T){
      abc <- coords(rocmodel, "b", ret = c("specificity", "sensitivity"))
      youdenres <- abc[1] + abc[2] - 1
      youdenres <- round(youdenres, digit)
      result <- c(group, subgroup1, other, youdenres)
      names(result) <- c("group","subgroup", retype, "youden")
    }else{
      result <- c(group, subgroup1,other)
      names(result) <- c("group", "subgroup",retype)
    }
  }
  return(result)
}

quiteROCFunc <- quietly(ROCStatFunc)

data("aSAH")
head(aSAH)
dat <- read.csv("dat", header = T)

### 计算outcome为结局变量的age的相关信息
quiteROCFunc(aSAH, group = "outcome", var = "age")$result

# Remove some NA values
dat <- dat[!is.na(dat$INDUCTION_RESPONSE),]

# 批量计算变量的ROC结果
## 定义group
Z10 <- read.csv("000B2_R_vs_NR_List.csv", header = T)

df <- read.table("000B_R_vs_NR_matrix2.txt", header = T)
df <- t(df)
df <- as.data.frame(df)
multigroup <- colnames(df)[1:ncol(df)]

df$ID <- row.names(df)
dat <- merge(Z10, df, by.x = "X", by.y = "ID")
###

dat <- dat[4:ncol(dat)]
rocRes <- lapply(multigroup, function(x) quiteROCFunc(dat, "Response", x)$result)

rocResDat <- do.call(rbind, rocRes)
rocResDat
rocResDat <- cbind(multigroup, rocResDat)
write.csv(rocResDat,"Gene_Only_ROC.csv")


###################
## https://www.rdocumentation.org/packages/pROC/versions/1.16.2/topics/plot.roc

Z9 <- read.csv("GSE16879_TOP1_Clinical.csv")
Z9 <- subset(Z0, Group.1 == "PRE.R" | Group.1 == "PRE.NR")

library(pROC)
X <- roc(Z6$GROUP1, Z6$V1287) # print the AUC (will contain the CI))

pdf("NUM1_ROC.pdf", 4.75,4.75)
plot.roc(A1_UP_DN$GROUP1 ,A1_UP_DN$V1606,
         levels=c("IFN.PRE.NR", "IFN.PRE.R"),
             legacy.axes=TRUE,
             ci=TRUE, 
             boot.n=10000,
             ci.alpha=0.95, 
             stratified=TRUE,
             auc.polygon=TRUE, 
             max.auc.polygon=FALSE, 
             grid=TRUE,
             print.auc=TRUE, 
             print.thres='best', 
             print.thres.best.method='y',
             print.thres.adj=c(-0.1, 1.25),
             print.thres.pattern="Cut-off: %.3f \n\nSp: %.3f \nSe: %.3f",
             print.thres.pattern.cex = 1,
             print.auc.cex = 1,
             print.auc.y=20,
             print.auc.x=80,
             cex.axis=1,
             cex.lab=1,
             print.thres.pch=16,
             print.thres.cex=1.25,
             cex.main=0.8,
         main="CD4.severe.asthma.HS.IVV.UP	+ \nGSE18893_CTRL_VS_TNF_TREATED_TREG_24H_DN", 
         percent=TRUE, 
         grid.col=c("green", "red"), 
         auc.polygon.col="lightblue", 
         col="black")
dev.off()

roc.menarche <- roc(formula = GROUP1 ~ V2522, 
                    levels=c("PRE.NR", "PRE.R"),
                    data = Z6)

roc.menarche

# bootstrap (delong)
ci(roc.menarche, of = "auc", method = "delong", boot.n=10000)    

ROC.Result <- ci.coords(roc.menarche, "best", 
               conf.level=0.95, 
               method = "DeLong", 
               boot.n=10000, 
               ret=c("threshold", "specificity", "sensitivity", "accuracy",
                     "tn", "tp", "fn", "fp", "npv", "ppv", "1-specificity",
                     "1-sensitivity", "1-accuracy", "1-npv", "1-ppv",
                     "precision", "recall"))
ROC.Result <- as.matrix(ROC.Result)