

setMethod("[", signature("factR"), function (x, i, j){
    # check features
    gtf <- methods::slot(x, "custom")
    genetxs <- methods::slot(x, name = "txdata")
    if(missing(i)){
        gtf <- gtf
    } else if(typeof(i) %in% c("integer", "double")){
        gtf <- gtf[i]
    } else if(typeof(i) %in% c("logical")){
        gtf <- gtf[i]
    } else {
        txs <- tryCatch(.getTxs(x, i),
                        error = function(e) rlang::abort("Feature not found"))
        gtf <- gtf[gtf$transcript_id %in% txs]
    }
    gtf[,j]
})
