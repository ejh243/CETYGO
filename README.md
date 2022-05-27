# CETYGO (CEll TYpe deconvolution GOodness)

The repository contains code to estimate cellular composition from bulk tissue using purified cell type reference DNA methylation data and to calculate a quality metric (DSRMSE) that reflects the accuracy of the this statistical deconvolution. CETYGO is applicable to any deconvolution algorithm and any tissue for which the relevant panel of reference data exists. In this initial application we have used Houseman's deconvolution algorithm and applied to mouse brain. We are currently working on additional applications for human whole blood and human brain.
For details on how to implement the method, see HowToApplyDSRMSE.R

Deconvolution models currently available are:
- Mouse (IN THE PROCESS OF VALIDATION)
    - Estimates proportion of cells that are NeuN+ve (neuronal) and NeuN-ve (non-neuronal).
	- Reference data generated using the custom Illumina array, HorvathMammalMethylChip40.
	
	
Scripts are adapted from minfi (http://bioconductor.org/packages/release/bioc/html/minfi.html) with the added error metric utility.


