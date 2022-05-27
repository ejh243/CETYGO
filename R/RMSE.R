#' Root mean square error (RMSE)
#'
#' Calculates the root mean squared error from two vectors of observed and expected data.
#' Assumes the two vectors are of the same length and that the pairs are the same order.
#' @param m A vector of expected values
#' @param o A vector sof observed values expected to be the same length as m

#' @return the root mean square error calculated from the two input vectors
#' @export

RMSE = function(m, o){
    sqrt(mean((m - o)^2))
}