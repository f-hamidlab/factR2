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

granges.factR <- function(object, ..., set = NULL, sort = TRUE){

    accepted.sets <- c("all", "gene", "transcript","exon",  "CDS",
                      "start_codon", "stop_codon", "AS")


    if(is.null(set)){
        set <- slot(object, "active.set")
    } else if(any(!set %in% accepted.sets)){
        rlang::warn("Incorrect set name or index, using active set")
        set <- slot(object, "active.set")
    }
    gtf <- methods::slot(object, "transcriptome")
    # get type of input

    if(!missing(...)){
        in.type <- S4Vectors::mcols(gtf)[,c("gene_id", "transcript_id", "gene_name")] %>%
            as.data.frame %>%
            tidyr::pivot_longer(gene_id:gene_name) %>%
            dplyr::filter(!is.na(value)) %>%
            dplyr::filter(value %in% c(...)) %>%
            dplyr::pull(name) %>%
            unique()
        gtf <- gtf[S4Vectors::mcols(gtf)[[in.type[1]]] %in% ...]
    }



    if("all" %in% set){
        return(gtf)
    } else {
        return(gtf[gtf$type %in% set])
    }
}
