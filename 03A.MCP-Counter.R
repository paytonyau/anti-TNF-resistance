#### xCELL and CIBERSORT were using the web tools

#### analysis on Deconvolution-To-Estimate-Immune-Cells followed by the authors' github insruction
#### https://github.com/holiday01/deconvolution-to-estimate-immune-cell-subsets

##### MCP-Counter #####
## https://github.com/ebecht/MCPcounter
## Install  R devtools and MCPcounter dependancy 'curl'
install.packages(c("devtools","curl")) 
library(devtools)
## Install the package
install_github("ebecht/MCPcounter",ref="master", subdir="Source")
library(MCPcounter)
## Save the matrix
expression = read.csv("merged_batch_effect_corrected_matrix", header = TRUE, row.names = 1)

#### example run
# MCPcounter.estimate(expression,
#                    featuresType=c("affy133P2_probesets","HUGO_symbols","ENTREZ_ID")[1],
#                    probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),
#                                               sep="\t",stringsAsFactors=FALSE,colClasses="character"),
#                    genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),
#                                           sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character", check.names=FALSE)
# )
#### end

## Run the program
Estimates = MCPcounter.estimate(expression, featuresType = "HUGO_symbols")
write.csv(Estimates, "MCPcounter_matrix.csv")

##### EPIC #####
### https://github.com/GfellerLab/EPIC
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
## Install the package
library(EPIC)

## Run the program
Estimates =  EPIC(expression)
## Save the matrix
write.csv(Estimates, "EPIC_matrix.csv")
