#' Perform simulations to test deconvolution model with constructed bulk
#' tissue profiles.
#'
#' Function to implement the Houseman deconvolution method for a provided
#'  set reference data and apply it to a series of  bulk profiles
#' constructed from user-provided proportions. The function first selects
#' the sites that will be used for the deconvolutions. It then generates
#' test bulk profiles which are weighted sums of reference profiles.
#' Finally it calculates estimates of the cellular composition and the
#' CETYGO score for each simulated bulk profile.
#'
#' We recommend that the training reference data and the test reference
#' data (which is used to construct bulk profiles for testing) are
#' distinct. Note that no normalisation is performed as part of this
#' function, it is recommended that data is thouroughly QC'd prior to
#' this analysis and your perfferred normalisation method applied.
#'
#' @param trainBetas A matrix of reference DNA methylation profiles,
#' where rows are sites and columns are samples. Note you need multiple
#' samples from the same cell type.
#' @param trainCellTypes A vector of which cell types to include in the
#' deconcolution model.
#' @param trainCellInd  A vector indicating which cell type each column in
#' trainBetas comes from.
#' @param testBetas A matrix with DNA methylation levels for reference cell
#' types to construct the bulk tissue profiles from. Format is one column per
#' cell type. Requires the same number of rows as trainBetas and in the same
#' order.
#' @param matrixSimProp A matrix of proportions to combine reference cell
#' types. Each row represents a different combination of cell types. Number
#' of columns must match the number of columns in testBetas, unless the last
#' column is the proportion of "Noise", and must be labelled as such.
#' @return A matrix with predicted cellular compositions and CETYGO score for
#' each simulated bulk tissue profile.
#' @export
#'
#' @examples
#' # create mean DNAm levels for 100 sites
#' set.seed(1327)
#' meanBetas <- runif(100, min = 0, max = 1)
#' # generate cell type diffs
#' meanCTDiff <- rnorm(100, mean = 0, sd = 0.2)
#' # create reference training data
#' refBetas <- cbind(
#'     matrix(meanBetas +
#'         rnorm(500, mean = 0, sd = 0.01),
#'     nrow = 100, byrow = FALSE
#'     ),
#'     matrix(meanBetas + rnorm(500, mean = 0, sd = 0.01) +
#'         meanCTDiff, nrow = 100, byrow = FALSE)
#' )
#' # force to lie between 0 and 1
#' refBetas[refBetas < 0] <- runif(sum(refBetas < 0), 0, 0.05)
#' refBetas[refBetas > 1] <- runif(sum(refBetas > 1), 0.95, 1)
#' rownames(refBetas) <- paste0("cg", 1:100)
#' # create test data
#' testBetas <- cbind(meanBetas, meanBetas + meanCTDiff) +
#'     rnorm(200, mean = 0, sd = 0.01)
#' # force to lie between 0 and 1
#' testBetas[testBetas < 0] <- runif(sum(testBetas < 0), 0, 0.05)
#' testBetas[testBetas > 1] <- runif(sum(testBetas > 1), 0.95, 1)
#' rownames(testBetas) <- paste0("cg", 1:100)
#' simProps <- matrix(data = c(0.5, 0.5, 0.2, 0.8), ncol = 2)
#' colnames(simProps) <- c("A", "B")
#'
#' runNoiseSimulations(
#'     trainBetas = refBetas,
#'     trainCellTypes = c("A", "B"),
#'     trainCellInd = c(rep("A", 5), rep("B", 5)),
#'     testBetas = testBetas,
#'     matrixSimProp = simProps
#' )
#'
runNoiseSimulations <- function(trainBetas, trainCellTypes,
                                trainCellInd, testBetas, matrixSimProp) {
    if (!identical(rownames(trainBetas), rownames(testBetas))) {
        stop("Rows of training and test data are not identical")
    }
    ## fit model
    model <- pickCompProbesMatrix(
        rawbetas = trainBetas,
        cellTypes = trainCellTypes,
        cellInd = trainCellInd,
        numProbes = 50,
        probeSelect = "auto"
    )

    ## create test data
    testBulkBetas <- createBulkProfiles(
        testBetas[rownames(model$coef), ],
        matrixSimProp
    )
    ## do cellular prediction with error
    predProp <- projectCellTypeWithError(testBulkBetas, model$coef)
    return(predProp)
}
