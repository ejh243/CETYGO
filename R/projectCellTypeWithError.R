#' Estimate cellular composition with error
#'
#' Estimate the cellular composition from a DNA methylation profile and
#' calculate the CETYGO score error metric associated with this estimate.
#'
#' @param YIN a matrix of DNA methylation levels from the samples that require
#' cell composition to be estimated.
#' @param coefs A matrix of deconvolution model parameters including
#' coefficients. This input can be generated with the
#' pickCompProbesMatrix()$coefs
#' @return A matrix with the estimated proportion of cell types, CETYGO and
#' the number of sites missing from the model.
#'
#' @examples
#' # using pre-trained model for whole blood provided with the package and
#' # sample whole blood profiles
#' rInd <- rownames(bulkdata)[rownames(bulkdata) %in% rownames(modelBloodCoef)]
#' predProp <- projectCellTypeWithError(bulkdata, modelBloodCoef[rInd, ])
#'
#' head(predProp)
#' @export
projectCellTypeWithError <- function(YIN, coefs) {

    ## if there's only 1 sample, pretend there are 2 so it doesn't break
    if (ncol(as.matrix(YIN)) == 1) {
        sampleDup <- 1
        YIN <- cbind(YIN, YIN)
    } else {
        sampleDup <- 0
    }
	
    ## subset the CpGs in YIN
	sharedProbes<-intersect(rownames(coefs), rownames(YIN))
	if(length(sharedProbes) < 5){
		stop(paste0("Only ", length(sharedProbes), " probes in common between trained model and bulk data. Stopping."))
		
	} else {
		if(length(sharedProbes) < 50){
		warning("Less than 50 probes in common between trained model and bulk data, this might affect the accuracy of the deconvolution")
		}
	}
    YIN <- YIN[sharedProbes, ]
    coefCellTypeIN <- coefs[sharedProbes,]


    projectCellType <- function(Y, coefCellType, contrastCellType = NULL,
                                nonnegative = TRUE, lessThanOne = FALSE) {
        if (is.null(contrastCellType)) {
            Xmat <- coefCellType
        } else {
            Xmat <- tcrossprod(coefCellType, contrastCellType)
        }

        nCol <- dim(Xmat)[2]
        if (nCol == 2) {
            obs <- which(apply(Y, 1, function(x) {
                return(sum(is.na(x)) == 0)
            }))
            Dmat <- crossprod(Xmat[obs, ])
            mixCoef <- t(apply(Y, 2, function(x) {
                solve(
                    Dmat,
                    crossprod(Xmat[obs, ], x[obs])
                )
            }))
            colnames(mixCoef) <- colnames(Xmat)
            return(mixCoef)
        } else {
            nSubj <- dim(Y)[2]

            mixCoef <- matrix(0, nSubj, nCol)
            rownames(mixCoef) <- colnames(Y)
            colnames(mixCoef) <- colnames(Xmat)

            if (nonnegative) {
                if (lessThanOne) {
                    Amat <- cbind(rep(-1, nCol), diag(nCol))
                    b0vec <- c(-1, rep(0, nCol))
                } else {
                    Amat <- diag(nCol)
                    b0vec <- rep(0, nCol)
                }
                for (i in seq_len(nSubj)) {
                    obs <- which(!is.na(Y[, i]))
                    Dmat <- crossprod(Xmat[obs, ])
                    mixCoef[i, ] <- solve.QP(
                        Dmat, crossprod(
                            Xmat[obs, ],
                            Y[obs, i]
                        ),
                        Amat, b0vec
                    )$sol
                }
            } else {
                for (i in seq_len(nSubj)) {
                    obs <- which(!is.na(Y[, i]))
                    Dmat <- crossprod(Xmat[obs, ])
                    mixCoef[i, ] <- solve(
                        Dmat,
                        t(Xmat[obs, ]) %*% Y[obs, i]
                    )
                }
            }
            return(mixCoef)
        }
    }

    ## error function
    getErrorPerSample <- function(applyIndex,
                                    predictedIN = mixCoef,
                                    coefDataIN = coefCellTypeIN,
                                    betasBulkIN = YIN) {
        trueBulk <- matrix(
            ncol = 1, nrow = nrow(coefDataIN),
            data = 0
        )
        for (i in seq_len(ncol(coefDataIN))) {
            trueBulk[, 1] <- trueBulk[, 1] +
                coefDataIN[, i] * predictedIN[applyIndex, i]
        }
        betasBulkIN <- t(apply(betasBulkIN, 1, function(x) {
            x[is.na(x)] <- 0
            return(x)
        }))
        error <- RMSE(trueBulk, betasBulkIN[, applyIndex])
        return(error)
    }

    mixCoef <- projectCellType(YIN, coefCellTypeIN)
    CETYGO <- sapply(seq_len(nrow(mixCoef)), getErrorPerSample)
	nMissingAll<-nrow(coefs) - nrow(coefCellTypeIN)
    nCGmissing <- apply(YIN, 2, function(x) {
        sum(is.na(x)) + nMissingAll
    })
    mixCoef <- cbind(mixCoef, CETYGO, nCGmissing)
    if (sampleDup == 1) {
        mixCoef <- mixCoef[, 1]
    }
    return(mixCoef)
}
