# CETYGO (CEll TYpe deconvolution GOodness)

The majority of epigenetic epidemiology studies to date have generated genome-wide profiles from bulk tissues (e.g. whole blood) however these are vulnerable to confounding from variation in cellular composition. Proxies for cellular composition can be mathematically derived from the bulk tissue profiles using a deconvolution algorithm however, there is no method to assess the validity of these estimates for a dataset where the true cellular proportions are unknown. We have developed an accuracy metric that quantifies the **CEll TYpe deconvolution GOodness (CETYGO)** score of a set of cellular heterogeneity variables derived from a genome-wide DNA methylation profile for an individual sample. The CETYGO score captures the deviation between a sample’s DNAm profile and its expected profile given the estimated cellular proportions and cell type reference profiles

This repository contains an R package with functions to estimate CETYGO. It builds on the existing functionality available through the minfi package, meaning it is straightforward to incorporate into existing pipelines and users can take advantage of the range of reference panels compatible with that methodology. 

## Citation

Vellame DS, Shireby G, MacCalman A, Dempster EL, Burrage J, Gorrie-Stone T, Schalkwyk LS, Mill J, Hannon E. Uncertainty quantification of reference-based cellular deconvolution algorithms. Epigenetics. 2023 Dec;18(1):2137659. doi: 10.1080/15592294.2022.2137659. Epub 2022 Dec 20. PMID: 36539387; PMCID: PMC9980651.

## Installation 

The [CETYGO package is available via GitHub](https://github.com/ds420/CETYGO) and can be installed using the devtools package. However, there are some pre-requistite packages that may need to be installed from Bioconductor first.

```
# install required packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("genefilter", "minfi"))
install.packages("quadprog")

# install devtools to install from GitHub
install.packages("devtools")
library(devtools)


install_github("ds420/CETYGO")
```

## Quick Start

Once you have successfully installed the CETYGO package you can recalculate your cell composition variables and associated CETYGO score. 

Within the CETYGO package we have provided functions and a a pre-trained model ```modelBloodCoef``` to 
enable the estimate of the composition of major blood cell types, as well as the CETYGO score
from a matrix of (normalised) beta values. We have also provided 10 
exemplar whole blood profiles generated with the 450K array in the R object 
```bulkdata```. Using these together we can quickly recalculate both the cellular 
proportions for six blood cell types and the CETYGO score for each sample. This 
can be done with the following code.

```

library(CETYGO)

rowIndex<-rownames(bulkdata)[rownames(bulkdata) %in% rownames(modelBloodCoef)]
predProp<-projectCellTypeWithError(bulkdata, modelBloodCoef[rowIndex,])

head(predProp)

```

For more details, and examples demonstrating how to implement the standard workflow including normalisation of the reference data with the bulk tissue (test) data, please see the vignette included with the package.  

## Interpretation of the CETYGO score

We have profiled the behaviour of the CETYGO score when applied to whole blood 
across a large number of empirical datasets and provide the following guidance 
for it's interpretation. 

1. CETYGO > 0.1 indicates the sample is not composed of the reference cell 
types, and is potentially the wrong tissue.

2. Elevated CETYGO can be indicative of a technically poor DNA methylation 
profile.

3. Purified cell types have higher CETYGO scores than bulk tissues.

4. Profiles generated with the EPIC array are associated with higher CETYGO 
scores than the 450K array.

5. Using our pre-trained model `modelBloodCoef`, across 3001 whole blood 
samples profiled with the 450K array, the median CETYGO score was 0.045  and 
the 95% "inter-quartile" range was 0.040-0.061. Across 3350 whole blood samples 
profiled with the EPIC array, the median CETYGO score was 0.057 and the  95% 
"inter-quartile" range was 0.050 - 0.069. We can use these results to propose 
an acceptable range of values. However, it is evident that this is technology 
and reference panel specific and therefore these boundaries may not be well 
calibrated for all applications. 


## Mathematical details

The CETYGO score captures the deviation between the observed DNAm profile and the expected profile for the given set of estimated cell type proportions. We can define the DNA methylation profile of a bulk tissue as the sum of DNA methylation levels measured in the constituent cell types weighted by the proportion of total cells represented by that cell type. 

**Equation 1**
$$B_{i,j}=\sum_{k=1}^N \left( p_{i,k} C_{i,j,k} \right)$$

where 

- $B_{i,j}$ represents the DNA methylation level in the bulk tissue for sample i at site j
- $p_{i,k}$ represents the proportion of cell type k in sample i
- $C_{i,j,k}$ represents the DNA methylation level for sample i at site j in cell type k, for N different cell types 
    
To estimate estimate $p_{i,k}$ for all (major) cell types, we can take advantage of reference profiles (i.e. $C_{i,j,k}$) available to the research community and methodologies, such as Houseman’s constraint projection approach, to solve for the unknown $p_{i,k}$. 

If the estimated cell proportions (denoted $\hat{p_{i,k}}$ ) are accurate then the expected bulk tissue profile given this composition of cell types should closely resemble the observed data. We can substitute our estimated cell proportions, ($\hat{p_{i,k}}$), back into Equation 1, to calculate the expected profile of DNA methylation values. 

**Equation 2**

$$\hat{B_{i,j}}=\sum_{k=1}^N \left( \hat{p_{i,k}} C_{i,j,k} \right)$$

CETYGO is then defined as the root mean square error (RMSE) between the observed bulk DNA methylation profile and the expected profile across the M cell type specific DNA methylation sites used to perform the deconvolution. 

**Equation 3**

$$ CETYGO_{i} = RMSE (B_{i}, \hat{B_{i}}) $$

## Brain reference panels

A tutorial on the brain reference panels that are available through CETYGO is [available here](https://github.com/ds420/CETYGO/wiki/Deconvolution-of-brain-cell-types)

## Additional package functionality

The CETYGO package also contains functions to generate constrcted bulk tissue profile from profiles of purified cell types at fixed proportions, both with and without noise. These may useful for testing out new reference panels. 

## Acknowledgements

We are grateful to the developers and contributors of the [minfi](https://github.com/hansenlab/minfi) package. By making their code open source and available through GitHub has enabled us to implement our metric within the existing framework  ultimately making it easier for users to add it to their workflow.
