#' Run factR workflow
#'
#' @description This wrapper function will perform the core factR workflow which
#' include (1) building of CDSs, (2) predicting of NMD-sensitive transcripts
#' and (3) labelling of alternative splicing segments.
#'
#' @param object factRObject
#' @param verbose whether to display workflow progress
#'
#' @return Updated factRObject
#' @export
#' @seealso \code{\link{buildCDS}}, \code{\link{predictNMD}}, \code{\link{findAltSplicing}}
#'
#' @rdname RunfactR
#' @examples
#' data(factRsample)
#' factRsample <- runfactR(factRsample)
setGeneric("runfactR", function(object, ...) standardGeneric("runfactR"))
setMethod("runfactR", "factR", function(
        object,
        NMD_threshold = 50,
        cons_db = "phastCons",
        verbose = TRUE) {

    object <- buildCDS(object, verbose)
    object <- predictNMD(object, NMD_threshold, verbose)
    object <- getAAsequence(object, verbose)
    object <- testASNMDevents(object, verbose)
    object <- getAScons(object, cons_db, "exon", 0)  # TODO: get upstream and downstream, with some padding
    object
})

















