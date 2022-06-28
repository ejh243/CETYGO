# CETYGO (CEll TYpe deconvolution GOodness)

 We have developed an accuracy metric that quantifies the **CEll TYpe deconvolution GOodness (CETYGO)** score of a set of cellular heterogeneity variables derived from a genome-wide DNA methylation profile for an individual sample. This R package provides users with functions to estimate these values, by building on the existing functionality available through the Biocpkg("minfi") package. While theorhetically the CETYGO score can be used in conjunction with any reference based deconvolution method, this package only contains code to calculate it in combination with Houseman's algorithm. By integrating our extension on top of the popular minfi package, means users can take advantage of the range of reference panels already formated for use, and add the calculation of the CETYGO score into their workflow with minimal changes to their existing scripts. 

## Citation

Uncertainty quantification of reference based cellular deconvolution algorithms
Dorothea Seiler Vellame, Gemma Shireby, Ailsa MacCalman, Emma L Dempster, Joe Burrage, Tyler Gorrie-Stone, Leonard S Schalkwyk, Jonathan Mill, Eilis Hannon
bioRxiv 2022.06.15.496235; doi: https://doi.org/10.1101/2022.06.15.496235

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

Once you have successfully installed the CETYGO package you can recalculate your cell composition variables and associated CETYGO score. Within the package we provide ```bulkdata``` an example set of whole blood DNA methylation profiles to demonstate how to use the package as follows. 

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

## Additional package functionality

The CETYGO package also contains functions to generate constrcted bulk tissue profile from profiles of purified cell types at fixed proportions, both with and without noise. These may useful for testing out new reference panels. 

## Acknowledgements

We are grateful to the developers and contributors of the [minfi](https://github.com/hansenlab/minfi) package. By making their code open source and available through GitHub has enabled us to implement our metric within the existing framework  ultimately making it easier for users to add it to their workflow.
