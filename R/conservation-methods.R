setGeneric("getAScons", function(
        object,
        db = "phastCons",
        type = "exon",
        padding = 50)  standardGeneric("getAScons"))
#' Quantify feature conservation
#'
#' @description Quantify sequence conservation scores of alternatively-spliced exons.
#'
#' @details
#' By default, this function quantifies the sequence conservation of the entire
#' exon together with 50 base-pair flanking sequence using the phastCons database.
#' Alternatively, this function can also quantify the mean conservation scores
#' of sequences flanking exons, or sequences upstream/downstream of alternative exons.
#'
#' The `padding` argument has different meaning depending on the `type` input.
#' If `type` is "exon", then `padding` refers to the amount of intronic sequences
#' on both sides of the exon to include in the quantification. If `type` is
#' "flanks", "upstream" or "downstream", then `padding` refers to the amount
#' of intronic sequnce to
#'
#' @param object factR object
#' @param db Database to extract sequence conservation. Can be "phastCons" (Default) or "phylop"
#' @param type Feature to quantify conservation. Can be one of the following:
#' \itemize{
#'  \item{"exon"}{ : Sequence conservation of the entire exon}
#'  \item{"intron"}{ : Sequence conservation of introns flanking exons. In the
#'  case of intron retention events, the ends of the sequence will be scored}
#'  \item{"flanks"}{ : Conservation of sequences flanking exons}
#'  \item{"upstream"}{ : Conservation of sequences upstream of exons}
#'  \item{"downstream"}{ : Conservation of sequences downstream of exons}
#' }
#' @param padding Additional width to pad the sequence by.
#' For cons_type "exons" and "flanks", padding will be added on both sides.
#'
#' @return factRObject with additional columns in AS metadata.
#'
#' The name of column(s) created will be "Cons.<type>.pad<padding>"
#'
#'
#' @export
#' @seealso \code{\link{runfactR}}
#'
#' @rdname getASCons
#' @examples
#' \dontrun{
#' ## Load sample factRObject
#' data(factRsample)
#'
#' ## By default, quantifies conservation of exon with 50bp paddings on each side
#' factRsample <- getAScons(factRsample)
#'
#' ## To quantify conservation of flanking intronic sequences
#' ### 100bp on both sides
#' factRsample <- getAScons(factRsample, type = "flanks", padding = 100)
#'
#' ### 50bp upstream
#' factRsample <- getAScons(factRsample, type = "upstream", padding = 50)
#'
#' ## To quantify conservation of sequences at the beginning and/or at the
#' end of an exon
#' ### 50bp Beginning and ending
#' factRsample <- getAScons(factRsample, type = "flanks", padding = -50)
#'
#' ### 50bp Beginningg
#' factRsample <- getAScons(factRsample, type = "upstream", padding = -50)
#'}
setMethod("getAScons", "factR", function(
          object,
          db = "phastCons",
          type = "exon",
          padding = 50) {

    return(.getASGScores(object, db, type, padding))
})

.getASGScores <- function(object, db = "phastCons", type = "exon", padding = 50){

    acceptedtypes <- c("exon","upstream", "downstream","flanks")
    if(!all(type %in% acceptedtypes)){
        unacceptedtype <- paste0(type[!type %in% acceptedtypes], collapse = ",")
        rlang::abort(stringr::str_glue("Input type(s) `{unacceptedtype}` not recognised"))
    }
    if((length(type) != length(padding)) & length(padding)>1){
        rlang::abort(stringr::str_glue(
            "Length of vector `padding` ({length(padding)}) not equal to `type` ({length(type)})"))
    }
    if(length(padding)==1){
        padding <- rep(padding, length(type))
    }


    .msgheader("Quantifiying conservation scores")
    consDB <- .GScorecheck(db, object@reference$build)
    if(is.null(consDB)){
        .msgsubwarn("Skipping conservation scoring (requirements not met)")
        return(object)
    }

    fullout <- do.call(cbind, lapply(seq_along(type), function(x){
        thistype <- type[x]
        thispadding <- padding[x]
        .msgsubinfo(stringr::str_glue(
            "Quantifying `{thistype}` conservation scores with {thispadding} padding"))
        AS.exons <- object@transcriptome[object@transcriptome$type %in% "AS"]
        AS.exons <- AS.exons %>% as.data.frame() %>%
            dplyr::distinct(AS_id, .keep_all = TRUE) %>%
            dplyr::filter(seqnames %in% GenomeInfoDb::seqnames(consDB)) %>%
            GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)



        ConsScoresDF <- tryCatch(
            {
                if(thistype == "flanks"){
                    AS.exons.start <- GenomicRanges::flank(AS.exons, start = TRUE, width = thispadding, both = F)
                    AS.exons.end <- GenomicRanges::flank(AS.exons, start = FALSE, width = thispadding, both = F)

                    AS.exons.start <- GenomicScores::gscores(consDB, AS.exons.start)
                    AS.exons.end <- GenomicScores::gscores(consDB, AS.exons.end)

                    scores <- rowMeans(cbind(AS.exons.start$default,
                                             AS.exons.end$default),
                                       na.rm = TRUE)
                    data.frame(scores = scores,
                               row.names = AS.exons$AS_id)
                } else if(thistype == "upstream"){
                    AS.exons.start <- GenomicRanges::flank(AS.exons, start = TRUE, width = thispadding, both = F)
                    AS.exons.start <- GenomicScores::gscores(consDB, AS.exons.start)

                    data.frame(scores = AS.exons.start$default,
                               row.names = AS.exons$AS_id)

                } else if(thistype == "downstream"){
                    AS.exons.end <- GenomicRanges::flank(AS.exons, start = FALSE, width = thispadding, both = F)
                    AS.exons.end <- GenomicScores::gscores(consDB, AS.exons.end)

                    data.frame(scores = AS.exons.end$default,
                               row.names = AS.exons$AS_id)

                } else if(thistype == "exon"){
                    scores <- GenomicScores::gscores(consDB, AS.exons+thispadding)
                    data.frame(scores = scores$default,
                               row.names = AS.exons$AS_id)
                } else {
                    rlang::abort()
                }
            },
            error = function(cond){
                rlang::abort("Unable to quantify exon conservation scores. Check if annotation package matches the genome used")
                return(AS.exons)
            }
        )
        ConsScoresDF$scores <- tidyr::replace_na(ConsScoresDF$scores, 0)

        # return df
        colname <- stringr::str_glue("Cons.{thistype}.pad{thispadding}")
        out <- data.frame(dat = ConsScoresDF$scores,
                          row.names = rownames(ConsScoresDF))
        colnames(out) <- colname
        return(out)
    }))


    # update ASE
    object <- addMeta(object, "AS", data = fullout)

    return(object)
}


.GScorecheck <- function(ConsScore, build){
    if("character" %in% is(ConsScore)){

        ## try getting GScore using listed genomes
        if(tolower(ConsScore) %in% c("phylop", "phastcons")){
            ConsScore <- tolower(ConsScore)
            if(build == "custom"){
                # raise warning and skip gscoring
                rlang::warn("Unable to determine GScore database from custom reference")
                return(NULL)
            } else{
                # TODO: check this upon adding other genomes
                data("genomes")
                db <- genomes[genomes$build == build,ConsScore][[1]]
                if(!is.na(db)){
                    .msgsubinfo(sprintf("Retrieving %s scores",
                                        db))
                    phast <- suppressMessages(GenomicScores::getGScores(db))
                    return(phast)
                } else {
                    # raise warning and skip gscoring
                    .msgsubwarn(stringr::str_glue("GScore database for {build} is unavailable"))
                    return(NULL)
                }
            }

        } else{
            ## try loading GScore package
            phast <- tryCatch(
                {
                    library(ConsScore, character.only = T )
                    .msgsubinfo(sprintf("Loaded %s package",
                                        ConsScore))
                    get(ConsScore)
                },
                error = function(cond){
                    GScoresList <- rownames(GenomicScores::availableGScores())
                    if(ConsScore %in% GScoresList){
                        .msgsubinfo(sprintf("Retrieving %s package",
                                            ConsScore))

                        suppressMessages(GenomicScores::getGScores(ConsScore))
                    } else {
                        rlang::abort(sprintf("%s score database not found",
                                             ConsScore))
                        return(NULL)
                    }
                })
            return(phast)
        }
    }
    else if("GScores" %in% is(ConsScore)){
        .msgsubinfo("Using loaded GScores object")
        return(ConsScore)
    } else {
        .msgsubwarn("`db` input not recognised")
        return(NULL)
    }
}
