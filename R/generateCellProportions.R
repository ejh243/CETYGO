
#' create matrix of cellular composition profiles consisting of proportion of noise and remaining compostion in the ratio of supplied cell types
#' 
#' @param cellProp A vector of the required cellular proportions. Must be named with cell labels.
#' @param noise A vector of the proportion of noise.
#' @return A matrix with length(noise) columns and length(cellProp)+1 columns.
#' @examples
#' bloodProp<-c(0.2,0.8)
#' names(bloodProp)<-c("A", "B")
#' generateCellProportions(bloodProp, c(0,0.5))
#' @export
generateCellProportions <- function(cellProp, noise){
	if(!is.null(names(cellProp))){
		cellWeight = cellProp/sum(cellProp) # standardarize  proportions to sum to 1
		matrixProp = (1 - noise) %o% cellWeight
		matrixProp = cbind(matrixProp, noise)
		colnames(matrixProp) = c(names(cellProp), "Noise")
		return(matrixProp)
	} else {
		stop("vector of cellular proportions does not have names")
	}
}