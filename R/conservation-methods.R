## TODO: do up documentation
### type can be exon, flanks, upstream, downstream

setGeneric("getAScons", function(object, ...) standardGeneric("getAScons"))
setMethod("getAScons", "factR", function(
          object,
          db = "phastCons",
          type = "exon",
          padding = 50) {

    return(.getASGScores(object, db, type, padding))
})


.getASGScores <- function(object, db = "phastCons", type = "exon", padding = 50){

    acceptedtypes <- c("exon","upstream", "downstream","flanks")
    if(!type %in% acceptedtypes){
        rlang::abort(stringr::str_glue("Input type `{type}` not recognised"))
    }

    .msgheader("Quantifiying conservation scores")
    consDB <- .GScorecheck(db, object@reference$build)
    if(is.null(consDB)){
        .msgsubwarn("Skipping conservation scoring (requirements not met)")
        return(object)
    }

    .msgsubinfo(stringr::str_glue(
        "Quantifying `{type}` conservation scores with {padding} padding"))
    AS.exons <- object@transcriptome[object@transcriptome$type %in% "AS"]
    AS.exons <- AS.exons %>% as.data.frame() %>%
        dplyr::distinct(AS_id, .keep_all = TRUE) %>%
        dplyr::filter(seqnames %in% GenomeInfoDb::seqnames(consDB)) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)



    ConsScoresDF <- tryCatch(
        {
            if(type == "flanks"){
                AS.exons.start <- GenomicRanges::flank(AS.exons, start = TRUE, width = padding, both = F)
                AS.exons.end <- GenomicRanges::flank(AS.exons, start = FALSE, width = padding, both = F)

                AS.exons.start <- GenomicScores::gscores(consDB, AS.exons.start)
                AS.exons.end <- GenomicScores::gscores(consDB, AS.exons.end)

                scores <- rowMeans(cbind(AS.exons.start$default,
                                                   AS.exons.end$default),
                                             na.rm = TRUE)
                data.frame(scores = scores,
                                  row.names = AS.exons$AS_id)
            } else if(type == "upstream"){
                AS.exons.start <- GenomicRanges::flank(AS.exons, start = TRUE, width = padding, both = F)
                AS.exons.start <- GenomicScores::gscores(consDB, AS.exons.start)

                data.frame(scores = AS.exons.start$default,
                           row.names = AS.exons$AS_id)

            } else if(type == "downstream"){
                AS.exons.end <- GenomicRanges::flank(AS.exons, start = FALSE, width = padding, both = F)
                AS.exons.end <- GenomicScores::gscores(consDB, AS.exons.end)

                data.frame(scores = AS.exons.end$default,
                           row.names = AS.exons$AS_id)

            } else if(type == "exon"){
                scores <- GenomicScores::gscores(consDB, AS.exons+padding)
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

    # update ASE
    colname <- stringr::str_glue("Cons.{type}.pad{padding}")
    object@sets$AS@rowData[colname] <- 0
    object@sets$AS@rowData[rownames(ConsScoresDF),colname] <- ConsScoresDF$scores

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
