#' Estimate cellular composition and associated error
#'
#' This is an adaptation of estimateCellCounts() from the minfi R package to
#' include the calculation of CETYGO alongside cellular deconvolution.
#' The arguments are as described in the minfi package. The estimation
#' of cellular composition is an implementaion of the Houseman et al
#' (2012) regression calibration approach algorithm to the Illumina 450k
#' microarray for deconvoluting heterogeneous tissue sources like blood.
#' For example, this function will take an RGChannelSet from a DNA
#' methylation (DNAm) study of blood, and return the relative proportions
#' of CD4+ and CD8+ T-cells, natural killer cells, monocytes, granulocytes,
#' and b-cells in each sample.
#'
#' The function currently supports cell composition estimation for blood,
#' cord blood, and the frontal cortex, through compositeCellType values of
#' "Blood", "CordBlood", and "DLPFC", respectively. Packages containing
#' the appropriate reference data should be installed before running the
#' function for the first time ("FlowSorted.Blood.450k",
#' "FlowSorted.DLPFC.450k", "FlowSorted.CordBlood.450k"). Each tissue
#' supports the estimation of different cell types, delimited via the
#' cellTypes argument. For blood, these are "Bcell", "CD4T", "CD8T",
#' "Eos", "Gran", "Mono", "Neu", and "NK" (though the default value for
#' cellTypes is often sufficient). For cord blood, these are "Bcell",
#' "CD4T", "CD8T", "Gran", "Mono", "Neu", and "nRBC". For frontal cortex,
#' these are "NeuN_neg" and "NeuN_pos". See documentation of individual
#' reference packages for more details.
#'
#' The meanPlot should be used to check for large batch effects in the data,
#' reducing the confidence placed in the composition estimates. This plot
#' depicts the average DNA methylation across the cell-type discrimating
#' probes in both the provided and sorted data. The means from the provided
#' heterogeneous samples should be within the range of the sorted samples.
#' If the sample means fall outside the range of the sorted means, the cell
#' type estimates will inflated to the closest cell type. Note that we
#' quantile normalize the sorted data with the provided data to reduce these
#' batch effects. DNA methylation of test samples,
#'
#' @param rgSet An input RGChannelSet  object with raw DNA methylation data
#' for the samples that require cell composition to be estimated.
#' @param compositeCellType Which composite cell type is being deconvoluted.
#' Should be one of "Blood", "CordBlood", or "DLPFC". See details.
#' @param processMethod How should the user and reference data be processed
#' together? Default input "auto" will use preprocessQuantile for Blood and
#' DLPFC and preprocessNoob otherwise, in line with the existing literature.
#' Set it to the name of a preprocessing function as a character if you want
#' to override it, like "preprocessFunnorm".
#' @param probeSelect How should probes be selected to distinguish cell types?
#' Options include "both", which selects an equal number (50) of probes (with
#' F-stat p-value < 1E-8) with the greatest magnitude of effect from the
#' hyper- and hypo-methylated sides, and "any", which selects the 100 probes
#' (with F-stat p-value < 1E-8) with the greatest magnitude of difference
#' regardless of direction of effect. Default input "auto" will use "any"
#' for cord blood and "both" otherwise, in line with previous versions of
#' this function and/or our recommendations. Please see the references for
#' more details.
#' @param cellTypes Which cell types, from the reference object, should be
#' we use for the deconvolution? See details.
#' @param referencePlatform The platform for the reference dataset; if the
#' input rgSet belongs to another platform, it will be converted using
#' convertArray.
#' @param returnAll Should the composition table and the normalized user
#' supplied data be return?
#' @param meanPlot Whether to plots the average DNA methylation across the
#' cell-type discrimating probes within the mixed and sorted samples.
#' @param verbose Should the function be verbose?
#' @param ... Passed to preprocessQuantile
#' @return A matrix with the estimated proportion of cell types across all
#' sites, CETYGO and number of sites missing from the model.
#' 
#' @export
estimateCellCountsWithError <- function(rgSet, compositeCellType = "Blood",
                                        processMethod = "auto", 
                                        probeSelect = "auto",
                                        cellTypes = c(
                                            "CD8T", "CD4T", "NK", "Bcell",
                                            "Mono", "Gran"
                                        ),
                                        referencePlatform = c(
                                            "IlluminaHumanMethylation450k",
                                            "IlluminaHumanMethylationEPIC",
                                            "IlluminaHumanMethylation27k"
                                        ),
                                        returnAll = FALSE, meanPlot = FALSE,
                                        verbose = TRUE, ...) {

    # Check inputs
    .isMatrixBackedOrStop(rgSet, "estimateCellCounts")
    .isRGOrStop(rgSet)
    rgSet <- as(rgSet, "RGChannelSet")
    referencePlatform <- match.arg(referencePlatform)
    rgPlatform <- sub(
        "IlluminaHumanMethylation",
        "",
        annotation(rgSet)[which(names(annotation(rgSet)) == "array")]
    )
    platform <- sub("IlluminaHumanMethylation", "", referencePlatform)
    if ((compositeCellType == "CordBlood") && (!"nRBC" %in% cellTypes)) {
        message("[estimateCellCountsWithError]
            Consider including 'nRBC' in argument 'cellTypes'
            for cord blood estimation.\n")
    }
    referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, platform)
    subverbose <- max(as.integer(verbose) - 1L, 0L)
    message(sprintf(
            "Loading reference data package for
            compositeCellType '%s' and
            referencePlatform '%s' (inferred package name is '%s').
            An error at this stage might mean the relevant package needs 
            to be installed",
            compositeCellType, platform, referencePkg))
    referenceRGset <- get(data(list = referencePkg))
    if (rgPlatform != platform) {
        rgSet <- convertArray(
            object = rgSet,
            outType = referencePlatform,
            verbose = subverbose
        )
    }
    if (!"CellType" %in% names(colData(referenceRGset))) {
        stop(
            sprintf("the reference sorted dataset (in this case '%s')
            needs to have a phenoData column called 'CellType'"),
            names(referencePkg)
        )
    }
    if (sum(colnames(rgSet) %in% colnames(referenceRGset)) > 0) {
        stop(
            "the sample/column names in the user set must not be in the ",
            "reference data "
        )
    }
    if (!all(cellTypes %in% referenceRGset$CellType)) {
        stop(sprintf(
            "all elements of argument 'cellTypes' needs to be part
            of the reference
            phenoData columns 'CellType' (containing the following
            elements: '%s')",
            paste(unique(referenceRGset$cellType), collapse = "', '")
        ))
    }
    if (length(unique(cellTypes)) < 2) {
        stop("At least 2 cell types must be provided.")
    }
    if ((processMethod == "auto") &&
        (compositeCellType %in% c("Blood", "DLPFC"))) {
        processMethod <- "preprocessQuantile"
    }
    if ((processMethod == "auto") &&
        (!compositeCellType %in% c("Blood", "DLPFC"))) {
        processMethod <- "preprocessNoob"
    }
    processMethod <- get(processMethod)
    if ((probeSelect == "auto") && (compositeCellType == "CordBlood")) {
        probeSelect <- "any"
    }
    if ((probeSelect == "auto") && (compositeCellType != "CordBlood")) {
        probeSelect <- "both"
    }

    if (verbose) {
        message(
            "[estimateCellCountsWithError]
                Combining user data with reference ",
            "(flow sorted) data.\n"
        )
    }
    newpd <- DataFrame(
        sampleNames = c(colnames(rgSet), colnames(referenceRGset)),
        studyIndex = rep(
            x = c("user", "reference"),
            times = c(ncol(rgSet), ncol(referenceRGset))
        ),
        stringsAsFactors = FALSE
    )
    referencePd <- colData(referenceRGset)
    combinedRGset <- combineArrays(
        object1 = rgSet,
        object2 = referenceRGset,
        outType = "IlluminaHumanMethylation450k"
    )
    colData(combinedRGset) <- newpd
    colnames(combinedRGset) <- newpd$sampleNames
    rm(referenceRGset)

    if (verbose) {
        message(
            "[estimateCellCountsWithError]
            Processing user and reference data ",
            "together.\n"
        )
    }
    if (compositeCellType == "CordBlood") {
        # NOTE: Here Shan wants to discard probes that they have decided
        #       shouldn't be used, for example multi-mapping probes. This is
        #       done by only using probes with names in the comptable.
        #       This is kind of ugly, and dataset dependent.
        combinedMset <- processMethod(combinedRGset, verbose = subverbose)
        compTable <- get(paste0(referencePkg, ".compTable"))
        combinedMset <- combinedMset[
            which(rownames(combinedMset) %in% rownames(compTable)),
        ]
    } else {
        combinedMset <- processMethod(combinedRGset)
    }
    rm(combinedRGset)

    # Extract normalized reference data
    referenceMset <- combinedMset[, combinedMset$studyIndex == "reference"]
    colData(referenceMset) <- as(referencePd, "DataFrame")
    mSet <- combinedMset[, combinedMset$studyIndex == "user"]
    colData(mSet) <- as(colData(rgSet), "DataFrame")
    rm(combinedMset)

    if (verbose) {
        message(
            "[estimateCellCountsWithError]
            Picking probes for composition ",
            "estimation.\n"
        )
    }
    compData <- pickCompProbes(
        mSet = referenceMset,
        cellTypes = cellTypes,
        compositeCellType = compositeCellType,
        probeSelect = probeSelect
    )
    coefs <- compData$coefEsts
    # TODO: Shouldn't be necessary to rm() anything
    rm(referenceMset)

    if (verbose) message("[estimateCellCountsWithError]
            Estimating composition and error.\n")
    counts <- projectCellTypeWithError(getBeta(mSet)[rownames(coefs), ], coefs)
    rownames(counts) <- colnames(rgSet)

    if (meanPlot) {
        smeans <- compData$sampleMeans
        smeans <- smeans[order(names(smeans))]
        sampleMeans <- c(
            colMeans2(
                x = getBeta(mSet),
                rows = match(rownames(coefs), rownames(mSet))
            ),
            smeans
        )
        sampleColors <- c(
            rep(1, ncol(mSet)),
            1 + as.numeric(factor(names(smeans)))
        )
        plot(sampleMeans, pch = 21, bg = sampleColors)
        legend("bottomleft",
            c("blood", levels(factor(names(smeans)))),
            col = seq_len(7),
            pch = 15
        )
    }
    if (returnAll) {
        return(list(
            counts = counts,
            compTable = compData$compTable,
            normalizedData = mSet
        ))
    } else {
        counts
    }
}

