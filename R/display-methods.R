#' @include factRObject-class.R

show.factR <- function(object){
    cat(sprintf("class: factRObject [version %s]\n", object@version))
    cat(sprintf("# transcriptome:\n   "))
    ngenes <- length(object[["gene"]]$gene_id)
    ntxs <- length(object[["transcript"]]$transcript_id)
    nnovel <- sum(object[["transcript"]]$novel == "yes")
    ncds <- sum(object[["transcript"]]$cds == "yes")
    cat(sprintf("%s genes; \n   ", ngenes))
    cat(sprintf("%s transcripts [%s novel]; \n   ", ntxs, nnovel))
    cat(sprintf("%s coding transcripts \n", ncds))
    cat(sprintf("# active set: %s\n", object@active.set))

    nsamples <- nrow(object@colData)
    samples <- rownames(object@colData)
    if(nsamples > 4){
        samples <- paste(c(samples[1], samples[2], " ... ", samples[nsamples-1],
                           samples[nsamples]), collapse = "  ")
    } else {
        samples <- paste(samples, collapse = ", ")
    }
    cat(sprintf("# samples (%s): %s\n", nsamples, samples))
}

granges.factR <- function(object, ..., set = NULL){
    if(is.null(set)){
        set <- slot(object, "active.set")
    } else if(!set %in% c(listSets(object), "all")){
        rlang::warn("Incorrect set name or index, using active set")
        set <- slot(object, "active.set")
    }
    gtf <- methods::slot(object, "transcriptome")
    out.type <- ifelse(set %in% c("transcript"), "transcript_id", "gene_id")
    feat <- .getFeat(object, ..., out = out.type)

    if(set == "all"){
        return(gtf[S4Vectors::mcols(gtf)[[out.type]] %in% feat])
    } else if(set == "transcript") {
        return(gtf[S4Vectors::mcols(gtf)[[out.type]] %in% feat & gtf$type %in% c("transcript", "exon")])
    } else {
        return(gtf[S4Vectors::mcols(gtf)[[out.type]] %in% feat & gtf$type %in% set])
    }
}
