.isMatrixBackedOrStop <- function(object, FUN) {
    if (!.isMatrixBacked(object)) {
        stop("'", FUN, "()' only supports matrix-backed minfi objects.",
             call. = FALSE)
    }
}

.isMatrixBacked <- function(object) {
    stopifnot(is(object, "SummarizedExperiment"))
    all(vapply(assays(object), is.matrix, logical(1L)))
}

.isRGOrStop <- function(object) {
    if (!is(object, "RGChannelSet")) {
        stop("object is of class '", class(object), "', but needs to be of ",
             "class 'RGChannelSet' or 'RGChannelSetExtended'")
    }
}
