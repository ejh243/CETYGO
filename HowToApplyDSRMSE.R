## Example script to use projectCellTypeWithError()

## available at github https://github.com/ds420/DSRMSE
## Dorothea Seiler Vellame ds420@exeter.ac.uk

## source the functions 
# replace PATHTOGITHUBFOLDERLOCALLY with own local path that github folder has been pulled to

setwd(PATHTOGITHUBFOLDERLOCALLY)
source("projectCellTypeWithError.R")

## to predict proportions in your own data use the following
projectCellType(YIN, model)

## YIN:
# should be a matrix of DNAm betas
# data should be class matrix 
# columns are samples
# rows are cgs with rownames in format cg********* 

## model: 
# should be one of the model options
# "mouseHumanHybrid" "mouse"

## Output:
# proportions for each cell type predicted
# error for prediction, calculated using DSRMSE
# number of CpGs missing (of the predictive CpGs)
