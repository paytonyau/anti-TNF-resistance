## https://github.com/ebecht/MCPcounter
##

install.packages(c("devtools","curl")) ##Installs devtools and the MCPcounter dependancy 'curl'
library(devtools)

install_github("ebecht/MCPcounter",ref="master", subdir="Source")

library(MCPcounter)

expression <- read.csv("000A_matrix.csv", header = T, row.names = 1)

## example
MCPcounter.estimate(expression,
                    featuresType=c("affy133P2_probesets","HUGO_symbols","ENTREZ_ID")[1],
                    probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character"),
                    genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",
                                     check.names=FALSE)
)
#### end

Estimates = MCPcounter.estimate(expression,featuresType="HUGO_symbols")


heatmap(as.matrix(ExampleEstimates),col=colorRampPalette(c("blue","white","red"))(100)) 
