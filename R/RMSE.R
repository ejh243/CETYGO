#' Root mean square error (RMSE)
#'
#' Calculates the root mean squared error from two vectors of observed and expected data.
#' Assumes the two vectors are of the same length and that the pairs are the same order.
#' @param m A vector of expected values
#' @param o A vector sof observed values expected to be the same length as m

#' @return the root mean square error calculated from the two input vectors
#' @export
#'
#' @examples
#' expected<-c(0.4,0.5,0.2)
#' observed<-c(0.38, 0.47, 0.25)

RMSE = function(m, o){
	if(length(m) == length(o)){
		return(sqrt(mean((m - o)^2)))
	} else {
		stop("Input vectors are not the same length")
	}
}