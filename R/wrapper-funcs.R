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
setGeneric("runfactR", function(object, verbose = TRUE) standardGeneric("runfactR"))
setMethod("runfactR", "factR", function(object, verbose = TRUE) {

    object <- buildCDS(object, verbose)
    object <- predictNMD(object, verbose = verbose)
    object <- getAAsequence(object, verbose)
    object <- testASNMDevents(object, verbose)
    object
})

















