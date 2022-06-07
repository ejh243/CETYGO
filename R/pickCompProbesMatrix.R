
#' Select probes for use in cellular deconvolution
#'
#' Using a DNA methylation dataset containing profiles for a panel of
#' reference cell types this function selects the sites for cellular
#' deconvolution and estimates the nessecary coefficients
#' This function is adapted from minfi pickCompProbes()
#' to take a matrix of beta values rather than mset
#'
#' @param rawbetas A matrix of DNA methylation levels from reference
#' cell types
#' @param cellInd A vector specifying which cell type each column of
#' rawbetas represents
#' @param cellTypes A vector of which cell types in cellInd to include
#' in generation of reference panel
#' @param numProbes The number of probes for each cell type to select
#' @param probeSelect How should probes be selected to distinguish
#' cell types? Options include
#' "both", which selects an equal number (50) of probes (with F-stat
#' p-value < 1E-8) with the greatest magnitude of effect from the
#' hyper- and hypo-methylated sides, and "any", which
#' selects the 100 probes (with F-stat p-value < 1E-8) with the
#' greatest magnitude of difference regardless of direction of effect.
#' @return A list with the estimate coefficients, summary of F tests
#' and cell type means.
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
#'     matrix(meanBetas + rnorm(500, mean = 0, sd = 0.01),
#'         nrow = 100, byrow = FALSE
#'     ),
#'     matrix(meanBetas + rnorm(500, mean = 0, sd = 0.01) + meanCTDiff,
#'         nrow = 100, byrow = FALSE
#'     )
#' )
#' # force to lie between 0 and 1
#' refBetas[refBetas < 0] <- runif(sum(refBetas < 0), 0, 0.05)
#' refBetas[refBetas > 1] <- runif(sum(refBetas > 1), 0.95, 1)
#' cellInd <- c(rep("A", 5), rep("B", 5))
#' rownames(refBetas) <- paste0("cg", seq_len(100))
#' pickCompProbesMatrix(
#'     rawbetas = refBetas,
#'     cellTypes = c("A", "B"),
#'     cellInd = cellInd, numProbes = 10
#' )
#'
pickCompProbesMatrix <- function(rawbetas, cellInd, cellTypes = NULL,
                            numProbes = 50, probeSelect = "both") {
    if (!is.null(cellTypes)) {
        if (!all(cellTypes %in% as.character(cellInd))) {
            stop("elements of argument 'cellTypes' is not part of 'cellInd'")
        }
        keep <- which(as.character(cellInd) %in% cellTypes)
        rawbetas <- rawbetas[, keep]
        cellInd <- cellInd[keep]
    }
    ## make cell type a factor
    cellInd <- factor(cellInd)
    ffComp <- rowFtests(rawbetas, cellInd)
    prof <- sapply(splitit(cellInd), function(i) rowMeans(rawbetas[, i]))
    r <- matrixStats::rowRanges(rawbetas)
    compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
    names(compTable)[1] <- "Fstat"
    names(compTable)[c(-2, -1, 0) + ncol(compTable)] <- c("low", "high", 
                                                           "range")
    tIndexes <- splitit(cellInd)
    tstatList <- lapply(tIndexes, function(i) {
        x <- rep(0, ncol(rawbetas))
        x[i] <- 1
        return(rowttests(rawbetas, factor(x)))
    })

    if (probeSelect == "any") {
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[, "p.value"] < 1e-8, ]
            yAny <- y[order(abs(y[, "dm"]), decreasing = TRUE), ]
            c(rownames(yAny)[seq_len(numProbes * 2)])
        })
    } else {
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[, "p.value"] < 1e-8, ]
            yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
            yDown <- y[order(y[, "dm"], decreasing = FALSE), ]
            c(rownames(yUp)[seq_len(numProbes)], 
               rownames(yDown)[seq_len(numProbes)])
        })
    }

    trainingProbes <- unique(na.omit(unlist(probeList)))
    if (length(grep("NA", trainingProbes)) > 0) {
        trainingProbes <- trainingProbes[-grep("NA", trainingProbes)]
    }
    rawbetas <- rawbetas[trainingProbes, ]

    pMeans <- colMeans(rawbetas)
    names(pMeans) <- cellInd

    form <- as.formula(sprintf(
        "y ~ %s - 1",
        paste(levels(cellInd), collapse = "+")
    ))
    phenoDF <- as.data.frame(model.matrix(~ cellInd - 1))
    colnames(phenoDF) <- sub("cellInd", "", colnames(phenoDF))
    if (ncol(phenoDF) == 2) { # two group solution
        X <- as.matrix(phenoDF)
        coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(rawbetas))
    } else { # > 2 group solution
        tmp <- validationCellType(
            Y = rawbetas, pheno = phenoDF,
            modelFix = form
        )
        coefEsts <- tmp$coefEsts
    }

    out <- list(
        coefEsts = coefEsts, compTable = compTable,
        sampleMeans = pMeans
    )
    return(out)
}

pickCompProbes <- function(mSet, cellTypes = NULL, numProbes = 50,
                    compositeCellType = compositeCellType,
                    probeSelect = probeSelect) {
    .isMatrixBackedOrStop(mSet)
    p <- getBeta(mSet)
    pd <- as.data.frame(colData(mSet))
    if (!is.null(cellTypes)) {
        if (!all(cellTypes %in% pd$CellType)) {
            stop(
                "elements of argument 'cellTypes' is not part of ",
                "'mSet$CellType'"
            )
        }
        keep <- which(pd$CellType %in% cellTypes)
        pd <- pd[keep, ]
        p <- p[, keep]
    }
    # NOTE: Make cell type a factor
    pd$CellType <- factor(pd$CellType, levels = cellTypes)
    ffComp <- rowFtests(p, pd$CellType)
    prof <- vapply(
        X = splitit(pd$CellType),
        FUN = function(j) rowMeans2(p, cols = j),
        FUN.VALUE = numeric(nrow(p))
    )
    r <- rowRanges(p)
    compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
    names(compTable)[1] <- "Fstat"
    names(compTable)[c(-2, -1, 0) + ncol(compTable)] <-
        c("low", "high", "range")
    tIndexes <- splitit(pd$CellType)
    tstatList <- lapply(tIndexes, function(i) {
        x <- rep(0, ncol(p))
        x[i] <- 1
        return(rowttests(p, factor(x)))
    })

    if (probeSelect == "any") {
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[, "p.value"] < 1e-8, ]
            yAny <- y[order(abs(y[, "dm"]), decreasing = TRUE), ]
            c(rownames(yAny)[seq(numProbes * 2)])
        })
    } else {
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[, "p.value"] < 1e-8, ]
            yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
            yDown <- y[order(y[, "dm"], decreasing = FALSE), ]
            c(
                rownames(yUp)[seq_len(numProbes)],
                rownames(yDown)[seq_len(numProbes)]
            )
        })
    }

    trainingProbes <- unique(unlist(probeList))
    p <- p[trainingProbes, ]

    pMeans <- colMeans2(p)
    names(pMeans) <- pd$CellType

    form <- as.formula(
        sprintf("y ~ %s - 1", paste(levels(pd$CellType), collapse = "+"))
    )
    phenoDF <- as.data.frame(model.matrix(~ pd$CellType - 1))
    colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
    if (ncol(phenoDF) == 2) {
        # Two group solution
        X <- as.matrix(phenoDF)
        coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
    } else {
        # > 2 groups solution
        tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
        coefEsts <- tmp$coefEsts
    }

    list(
        coefEsts = coefEsts,
        compTable = compTable,
        sampleMeans = pMeans
    )
}
