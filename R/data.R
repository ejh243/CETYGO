#' Whole blood DNA methylation profiles.
#'
#' A sample dataset of whole blood DNA methylation profiles from 10 individuals
#'
#' @format A matrix with 599 rows and 10 variables
"bulkdata"

#' Cortical DNA methylation profiles.
#'
#' A sample dataset of cortical DNA methylation profiles from 10 individuals
#'
#' @format A matrix with 1862 rows and 10 variables
"pfcdata"

#' Coefficients for 6 blood cell types
#'
#' Coefficients for blood cell types for use in deconvoluting cell composition
#' from whole blood profiles. These data were generated using the
#' pickCompProbes() function from the minfi package and DNA methylation data
#' from the FlowSorted.Blood.450k R package.
#'
#' @format A matrix with 600 rows and 6 variables
"modelBloodCoef"

#' Coefficients for brain cell types
#'
#' Coefficients for brain cell types for use in deconvoluting cell composition
#' from bulk brain profiles. These data were generated using the either the
#' pickCompProbes() function from the minfi package or the IDOLoptimize() 
#' function from the IDOL package and DNA methylation data
#' from GSE234520. Each matrix represents a different
#' panel of cell types
#'
#' @format A list with 2 lists, where each element is a further list of matrices
"modelBrainCoef"


