#' Generate bulk tissue profiles from purified cell types
#'
#' createBulkProfiles
#' @param purifiedBetas 
#' @param matrixSimProp 
#' @return 
#' @export
createBulkProfiles<-function(purifiedBetas, matrixSimProp){
	## if the last column is entitled noise, add a random noise vector
	if(colnames(matrixSimProp)[ncol(matrixSimProp)] == "Noise"){
		## add noise from uniform distribution
		purifiedBetas <- cbind(purifiedBetas, runif(nrow(purifiedBetas), min = 0, max = 1))
	}
	testBulkBetas <- purifiedBetas %*% t(matrixSimProp)
	return(testBulkBetas)
}