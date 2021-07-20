## Example script to use projectCellTypeWithError()

## available at github https://github.com/ds420/DSRMSE
## Dorothea Seiler Vellame ds420@exeter.ac.uk

## source the functions 
# replace PathToGithubFolderLocally with own local path that github folder has been pulled to

setwd(PathToGithubFolderLocally)
source("projectCellTypeWithError.R")

## to predict proportions in your own data use the following
projectCellType(YIN, modelType = "ownModel", ownModelData = model)

## YIN:
# should be a matrix of DNAm betas
# data should be class matrix 
# columns are samples
# rows are cgs with rownames in format cg********* 

## model: 
# should be one of the model options
# "mouseHumanHybrid" "mouse"
# if modelType "ownModel" is used, the ownModelData should be in the same format as the other model inputs
# a list with the first list object containing a matrix of CpGs and model coeficients per cell type
# model can be created using pickCompProbes

## Output:
# proportions for each cell type predicted
# error for prediction, calculated using DSRMSE
# number of CpGs missing (of the predictive CpGs)


## To create a deconvolution model from your own reference data, use
pickCompProbes(rawbetas, cellInd, cellTypes = NULL, numProbes = 50, probeSelect = probeSelectIN)

## where rawbetas are the betas matrix: 
# columns are samples
# rows are cpgs with rownames in format cg********* 
# cellInd: unique names of cell types available 
# cellTypes: a factor list of cell types
# numProbes: the number of cpgs used toe predict per cell type

