# Deconvolution models with Deconvolution Specific Root Mean Squared Error (DSRMSE)

For usage of models, see HowToApplyDSRMSE.R

Deconvolution models currently available are:
- Mouse 
    - Created in custom Illumina array, HorvathMammalMethylChip40. 
    - Created using prefrontal cortex.
    - Differentiates between neuronal (NeuN+) and non-neuronal (NeuN-) cell types.
- Mouse human hybrid 
    - Created in custom Illumina array, HorvathMammalMethylChip40, also applicable to Illumina EPIC (850k) array. 
    - Created using prefrontal cortex.
    - Differentiates between neuronal (NeuN+) and non-neuronal (NeuN-) cell types.
- Own model
    - You can apply your own reference data, which should be in the same format as the above models containing cpg coefficients, containing cpg id's as rownames, cell types as columns.
    - The deconvolution model can be made using function pickCompProbes.R

Coming soon:
- Human blood 
  - Created in Illumina 450k, applicable to 450k and EPIC arrays.
  - Differentiates between CD8 T cells, CD4 T cells, natural killer cells,  B cells, monocytes, and neutrophils.
- Human brain
  - Created in post mortem prefrontal cortex, profiled in EPIC, applicable to 450k and EPIC arrays.
  - Differentiates between neurons (NeuN+), oligodendrocytes (Sox10+), glia (IRF8+) and other cell types.
