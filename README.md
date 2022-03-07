# PDO_BME_EKGel

## Introduction
This is a project to regenerate the RNAseq compariosn between PDOs clutured in two mediums  in this paper:
"Biomimetic hydrogel supports initiation and growth
of patient-derived breast tumor organoids"
citation's info to be added

The script for the analysis is written in R.


----

## Dependencies:
These R packages need to be installed in order to run the analysis script:
- Biobase
- ggrepel
- ggplot2
- ComplexHeatmap
- circlize
- DESeq2
- snowfall
- GSA
- piano


----
## Reproducibility of the Analysis:
- Once the project is downloaded to the user computer, the user needs to navigate to the main directory of the project "PDO_BME_EKGel-master".
- Inside the main directory, there is an R script file named "analysis.R". Running this script will regenerate the RNAseq compariosn between PDOs clutured in two mediums  in this paper.

**Important Note:** the user needs to set the working directory inside the script file before running it, i.e. change the following code inside the script:

`setwd("PDO_BME_EKGel") # set to path of cloned project in your machine`

- Included in the project are:
1- RNAseq processed data for all PDOs clutured in different mediums, `Data/UHN_BR_PDOs_hg38_kallisto_mixed.rda`.
2- REACTOME genesets used for the pathway analysis, `Data/c2.cp.reactome.v7.4.symbols.gmt`.
