# CETYGO (CEll TYpe deconvolution GOodness)

 We have developed an accuracy metric that quantifies the **CEll TYpe deconvolution GOodness (CETYGO)** score of a set of cellular heterogeneity variables derived from a genome-wide DNA methylation profile for an individual sample. This R package provides users with functions to estimate these values, by building on the existing functionality available through the Biocpkg("minfi") package. While theorhetically the CETYGO score can be used in conjunction with any reference based deconvolution method, this package only contains code to calculate it in combination with Houseman's algorithm. By integrating our extension on top of the popular minfi package, means users can take advantage of the range of reference panels already formated for use, and add the calculation of the CETYGO score into there workflow with minimal input. 

## Citation

## Installation 

The [CETYGO package is available via GitHub](https://github.com/ds420/CETYGO) and can be installed using the devtools package. However, there are some pre-requistite packages that need to be installed from Bioconductor first.

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


## Additional package functionality

The CETYGO package also contains functions to generate constrcted bulk tissue profile from profiles of purified cell types at fixed proportions, both with and without noise. These may useful for testing out new reference panels. 

## Acknowledgements

We are grateful to the developers and contributors of the [minfi](https://github.com/hansenlab/minfi) package. By making their code open source and available through GitHub has enabled us to implement our metric within the existing framework  ultimately making it easier for users to add it to their workflow.
