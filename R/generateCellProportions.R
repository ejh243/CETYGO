
#' Create matrix of cellular composition profiles consisting of proportion
#'  of noise and remaining compostion in the ratio of supplied cell types
#'
#' @param cellProp A vector of the required cellular proportions. Must be
#' named with cell labels.
#' @param noise A vector of the proportion of noise.
#' @return A matrix with length(noise) columns and length(cellProp)+1
#' columns.
#' @export
#' @examples
#' meanBloodProp <- c(3.01, 13.4, 6.13)
#' names(meanBloodProp) <- c("Bcell", "CDT4+", "CDT8")
#' noise <- c(seq(0, 0.1, 0.02))
#' generateCellProportions(meanBloodProp, noise)
#'
generateCellProportions <- function(cellProp, noise) {
    if (!is.null(names(cellProp))) {
        # standardarize  proportions to sum to 1
        cellWeight <- cellProp / sum(cellProp)
        matrixProp <- (1 - noise) %o% cellWeight
        matrixProp <- cbind(matrixProp, noise)
        colnames(matrixProp) <- c(names(cellProp), "Noise")
        return(matrixProp)
    } else {
        stop("vector of cellular proportions does not have names")
    }
}
