#' Run factR workflow
#'
#' @description This wrapper function will perform the core factR workflow which
#' includes:
#' \itemize{
#'  \item{1: }{ Building CDS structures }
#'  \item{2: }{ Predicting NMD-sensitive transcripts }
#'  \item{3: }{ Translating open-reading frames }
#'  \item{4: }{ Quantify sequence conservation of alternative exons }
#' }
#'
#' @param object factR object
#' @param NMD_threshold Minimum distance between PTC and downstream exon-exon junction to trigger NMD (Default: 50)
#' @param cons_db Database to extract sequence conservation. Can be "phastCons" (Default) or "phylop"
#' @param cons_type Feature to quantify conservation. Can be one of the following:
#' \itemize{
#'  \item{"exon"}{ : Sequence conservation of the entire exon}
#'  \item{"flanks"}{ : Conservation of sequences flanking exons}
#'  \item{"upstream"}{ : Conservation of sequences upstream of exons}
#'  \item{"downstream"}{ : Conservation of sequences downstream of exons}
#' }
#' @param cons_padding Additional width to pad the sequence by. For cons_type "exons" and "flanks", paddings will be added on both sides.
#' @param verbose Whether to show pipeline messages (Default: TRUE)
#'
#' @return Update factRObject with additional data from runfactR
#' @export
#' @seealso \code{\link{buildCDS}}, \code{\link{predictNMD}}, \code{\link{testASNMDevents}}, \code{\link{getAScons}}, \code{\link{getAAsequence}}
#'
#' @rdname RunfactR
#' @examples
#' data(factRsample)
#' factRsample <- runfactR(factRsample)
setGeneric("runfactR", function(
        object,
        NMD_threshold = 50,
        cons_db = "phastCons",
        cons_type = "exon",
        cons_padding = 0,
        verbose = TRUE) standardGeneric("runfactR"))

setMethod("runfactR", "factR", function(
        object,
        NMD_threshold = 50,
        cons_db = "phastCons",
        cons_type = "exon",
        cons_padding = 0,
        verbose = TRUE) {

    object <- buildCDS(object, verbose)
    object <- predictNMD(object, NMD_threshold, verbose)
    object <- getAAsequence(object, verbose)
    object <- testASNMDevents(object, verbose)
    object <- getAScons(object, cons_db, cons_type, cons_padding)  # TODO: get upstream and downstream, with some padding
    object
})

















