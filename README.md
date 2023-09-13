
### Hyperactive neutrophil chemotaxis contributes to anti-tumor necrosis factor-α treatment resistance in inflammatory bowel disease

Anti-TNF-α treatment in inflammatory bowel disease has been the subject of extensive research, with the objective of elucidating molecular pathways and establishing reliable diagnostic biomarkers for patient response. The investigation revealed no significant variances in immune microenvironment scores between responders and non-responders to anti-TNF-α treatment. However, an augmented presence of neutrophils, endothelial cells, and B-cells was observed in non-responders at baseline, suggesting that chemotaxis pathways may play a pivotal role in mediating resistance to treatment. Notably, the study identified Interleukin 13 receptor subunit alpha 2 (IL13RA2) as a potential biomarker with a promising capacity to predict anti-TNF-α treatment response, exhibiting a sensitivity of 68.13% and specificity of 84.93%. Consequently, the study concluded that hyperactive neutrophil chemotaxis could potentially modulate responses to anti-TNF-α treatment.

#### Table of contents
##### A. IBD patients data with anti-TNFα treatment response ([GSE16879](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16879), [GSE23597](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23597), [GSE52746](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52746), [GSE73661](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73661) & [GSE92415](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92415))

[01. Expression microarray data normalisation](https://github.com/paytonyau/anti-TNF-resistance/blob/main/01.microarray.normalisation.R)

[02. Data integration and batch effect correction](https://github.com/paytonyau/anti-TNF-resistance/blob/main/02.ComBat.SVA.R)

[03. In silico flow cytometry analysis](https://github.com/paytonyau/anti-TNF-resistance/blob/main/03.cellular_deconvolution.R)

[04. t-test Loop](https://github.com/paytonyau/anti-TNF-resistance/blob/main/04_T-Test.R)

[05. Three-dimensional plot](https://github.com/paytonyau/anti-TNF-resistance/blob/main/05_3D_Plot.R)

[06. Differentially expressed genes analysis](https://github.com/paytonyau/anti-TNF-resistance/blob/main/06_DEGs.R)

[07. ROC loop](https://github.com/paytonyau/anti-TNF-resistance/blob/main/07_ROC.curves.R)

[*Normalised matrix](https://github.com/paytonyau/anti-TNF-resistance/blob/main/IBD_matrix.7z)

[*MetaScape Result](https://github.com/paytonyau/anti-TNF-resistance/blob/main/IBD_MetaScape_outputs.zip)

##### B. Expermintal data on lipopolysaccharide (LPS)-induced Neutrophils ([GSE145918](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145918))

[08. RNA-seq data normalisation](https://github.com/paytonyau/anti-TNF-resistance/blob/main/08_LPS_Neutrophils_RNA-seq_Norm.R)

[09. Matched GO term calculation](https://github.com/paytonyau/anti-TNF-resistance/blob/main/09_LPS_Neutrophils_Groups_Comparsions.R)


Detailed explanation can be found in the manuscript below

Yau, Tung On, et al. [Hyperactive neutrophil chemotaxis contributes to anti‐tumor necrosis factor‐α treatment resistance in inflammatory bowel disease.](https://doi.org/10.1111/jgh.15764) Journal of Gastroenterology and Hepatology 37.3 (2022): 531-541.
