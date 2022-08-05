.getASGScores <- function(object, db, type = "flanks", padding = 200){

    # TODO: check for input `type`

    consDB <- .GScorecheck(db)

    rlang::inform("Quantifying exon conservation scores")
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
            } else{
                scores <- GenomicScores::gscores(consDB, AS.exons+padding)
                data.frame(scores = scores$default,
                                  row.names = AS.exons$AS_id)
            }
        },
        error = function(cond){
            rlang::abort("Unable to quantify exon conservation scores. Check if annotation package matches the genome used")
            return(AS.exons)
        }
    )

    # update granges
    object@transcriptome$score <- ifelse(!is.na(object@transcriptome$AS_id),
                                         ConsScoresDF[object@transcriptome$AS_id,]$scores,
                                         NA)

}


.GScorecheck <- function(ConsScore){
    if("character" %in% is(ConsScore)){

        ## try loading GScore package
        phast <- tryCatch(
            {
                library(ConsScore, character.only = T )
                rlang::inform(sprintf("Loaded %s package",
                                      ConsScore))
                get(ConsScore)
            },
            error = function(cond){
                GScoresList <- rownames(GenomicScores::availableGScores())
                if(ConsScore %in% GScoresList){
                    rlang::inform(sprintf("Retrieving %s scores",
                                          ConsScore))
                    GenomicScores::getGScores(ConsScore)
                } else {
                    rlang::abort(sprintf("%s score database not found",
                                          ConsScore))
                    return(NULL)
                }
            })
        return(phast)
    }
    else if("GScores" %in% is(ConsScore)){
        rlang::inform("Using loaded GScores object")
        return(ConsScore)
    } else {
        rlang::abort("`db` input not recognised")
    }

}
