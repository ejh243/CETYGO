.isMatrixBackedOrStop <- function(object, FUN) {
    if (!.isMatrixBacked(object)) {
        stop("'", FUN, "()' only supports matrix-backed minfi objects.",
            call. = FALSE
        )
    }
}

.isMatrixBacked <- function(object) {
    stopifnot(is(object, "SummarizedExperiment"))
    all(vapply(assays(object), is.matrix, logical(1L)))
}

.isRGOrStop <- function(object) {
    if (!is(object, "RGChannelSet")) {
        stop(
            "object is of class '", class(object),
            "', but needs to be of ",
            "class 'RGChannelSet' or 'RGChannelSetExtended'"
        )
    }
}

validationCellType <- function(Y, pheno, modelFix, modelBatch = NULL,
                                L.forFstat = NULL, verbose = FALSE) {
    N <- dim(pheno)[1]
    pheno$y <- rep(0, N)
    xTest <- model.matrix(modelFix, pheno)
    sizeModel <- dim(xTest)[2]
    M <- dim(Y)[1]
    
    if (is.null(L.forFstat)) {
        # NOTE: All non-intercept coefficients
        L.forFstat <- diag(sizeModel)[-1, ]
        colnames(L.forFstat) <- colnames(xTest)
        rownames(L.forFstat) <- colnames(xTest)[-1]
    }
    
    # Initialize various containers
    sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
    coefEsts <- matrix(NA, M, sizeModel)
    coefVcovs <- list()
    
    if (verbose) message("[validationCellType] ")
    # Loop over each CpG
    for (j in seq_len(M)) {
        # Remove missing methylation values
        ii <- !is.na(Y[j, ])
        nObserved[j] <- sum(ii)
        pheno$y <- Y[j, ]
        
        if (j %% round(M / 10) == 0 && verbose) message(".") # Report progress
        
        # Try to fit a mixed model to adjust for plate
        try({
            if (!is.null(modelBatch)) {
                fit <- try(
                    lme(modelFix, random = modelBatch, data = pheno[ii, ])
                )
                # NOTE: If LME can't be fit, just use OLS
                OLS <- inherits(fit, "try-error")
            } else {
                OLS <- TRUE
            }
            
            if (OLS) {
                fit <- lm(modelFix, data = pheno[ii, ])
                fitCoef <- fit$coef
                sigmaResid[j] <- summary(fit)$sigma
                sigmaIcept[j] <- 0
                nClusters[j] <- 0
            } else {
                fitCoef <- fit$coef$fixed
                sigmaResid[j] <- fit$sigma
                sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
                nClusters[j] <- length(fit$coef$random[[1]])
            }
            coefEsts[j, ] <- fitCoef
            coefVcovs[[j]] <- vcov(fit)
            
            useCoef <- L.forFstat %*% fitCoef
            useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
            Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef)) / sizeModel
        })
    }
    if (verbose) message(" done\n")
    
    # Name the rows so that they can be easily matched to the target data set
    rownames(coefEsts) <- rownames(Y)
    colnames(coefEsts) <- names(fitCoef)
    degFree <- nObserved - nClusters - sizeModel + 1
    
    # Get P values corresponding to F statistics
    Pval <- 1 - pf(Fstat, sizeModel, degFree)
    
    list(
        coefEsts = coefEsts,
        coefVcovs = coefVcovs,
        modelFix = modelFix,
        modelBatch = modelBatch,
        sigmaIcept = sigmaIcept,
        sigmaResid = sigmaResid,
        L.forFstat = L.forFstat,
        Pval = Pval,
        orderFstat = order(-Fstat),
        Fstat = Fstat,
        nClusters = nClusters,
        nObserved = nObserved,
        degFree = degFree
    )
}

splitit <- function(x) {
    split(seq(along = x), x)
}
