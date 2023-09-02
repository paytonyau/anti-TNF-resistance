#### (0) ESTIMATE (Estimation of STromal and Immune cells in MAlignant Tumor tissues using Expression data) https://bioinformatics.mdanderson.org/public-software/estimate/

#### (1) xCELL and (2) CIBERSORT were using the web available tools

#### (3) analysis on Deconvolution-To-Estimate-Immune-Cells followed by the authors' github insruction
#### https://github.com/holiday01/deconvolution-to-estimate-immune-cell-subsets

#### (4) MCP-Counter ####
## https://github.com/ebecht/MCPcounter
## Install  R devtools and MCPcounter dependancy 'curl'
install.packages(c("devtools","curl")) 
library(devtools)

## Install the package
install_github("ebecht/MCPcounter",ref="master", subdir="Source")
library(MCPcounter)

## Save the matrix
expression = read.csv("merged_batch_effect_corrected_matrix", header = TRUE, row.names = 1)

## Run the program
Estimates = MCPcounter.estimate(expression, featuresType = "HUGO_symbols")

# Save the matrix
write.csv(Estimates, "MCPcounter_matrix.csv")

#### (5) EPIC #####
### https://github.com/GfellerLab/EPIC
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)

# Install the package
library(EPIC)

# Run the program
Estimates =  EPIC(expression)

# Save the matrix
write.csv(Estimates, "EPIC_matrix.csv")