#' Generate bulk tissue profiles from purified cell types
#'
#' Given some reference DNA methylation profiles, this function combines them
#' in a weighted sum using the user-provided proportions to construct bulk
#' tissue profiles in known ratios. There is also the option to specifiy the
#' addition of a proportion of noise. If this is required a random vector of
#' DNA methylation levels is generated from a uniform distribution.
#'
#'
#' @param purifiedBetas A matrix with DNA methylation levels for reference
#' cell types. Format is one column per cell type.
#' @param matrixSimProp A matrix of proportions to weight reference cell
#' types. Each row represents a different combination of cell types.
#' Number of columns must match the number of columns in purifiedBetas,
#' unless the last column is the proportion of noise, and must be labelled
#'  as such. Columns must be named.
#'
#' @return A matrix of DNA methylation levels, with the same number of rows
#' as rows in purifiedBetas and the same number of columns as rows in
#' matrixSimProp.
#'
#' @examples
#' # 3 cell type example
#' cellProps <- matrix(c(0.1, 0.4, 0.5, 0.2, 0.2, 0.6), ncol = 3)
#' colnames(cellProps) <- c("A", "B", "C")
#' refBetas <- matrix(runif(10 * 3, min = 0, max = 1), ncol = 3)
#' createBulkProfiles(refBetas, cellProps)
#'
#' # add noise to example
#' cellProps <- cbind(cellProps, c(0.1, 0.2))
#' colnames(cellProps)[4] <- "Noise"
#' createBulkProfiles(refBetas, cellProps)
#' @export
createBulkProfiles <- function(purifiedBetas, matrixSimProp) {
    ## if the last column is entitled noise, add a random noise vector
    if (is.null(colnames(matrixSimProp))) {
        stop("Column names of cell type proportions must be named")
    }
    if (colnames(matrixSimProp)[ncol(matrixSimProp)] == "Noise") {
        ## add noise from uniform distribution
        purifiedBetas <- cbind(
            purifiedBetas,
            runif(nrow(purifiedBetas), min = 0, max = 1)
        )
    }
    if (ncol(purifiedBetas) == ncol(matrixSimProp)) {
        testBulkBetas <- purifiedBetas %*% t(matrixSimProp)
        return(testBulkBetas)
    } else {
        stop("Number of cell types in input matrix do not match")
    }
}
