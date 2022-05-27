#' Perform simulations to test deconvolution with proportion of noise
#'
#' createBulkProfiles
#' @param trainBetas 
#' @param trainCellTypes 
#' @param trainCellInd 
#' @param testBetas 
#' @param matrixSimProp 
#' @return 
#' @export

runNoiseSimulations<-function(trainBetas, trainCellTypes, trainCellInd, testBetas, matrixSimProp){
	## fit model
	model <- pickCompProbesMatrix(rawbetas = trainBetas,
									   cellTypes = trainCellTypes,
									   cellInd = trainCellInd,
									   numProbes = 50,
									   probeSelect = "auto")

	## create test data
	testBulkBetas <-createBulkProfiles(testBetas[rownames(model$coef),], matrixSimProp)
	## do cellular prediction with error
	predProp<-projectCellTypeWithError(testBulkBetas, model$coef)
	return(predProp)
}

