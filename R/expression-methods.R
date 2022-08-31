
setGeneric("addTxCounts", function(object, countData, sampleData,
                                   verbose = TRUE) standardGeneric("addTxCounts"))
setMethod("addTxCounts", "factR", function(object, countData, sampleData, verbose = TRUE) {

    # catch missing args
    mandargs <- c("countData", "sampleData")
    passed <- names(as.list(match.call())[-1])
    if (any(!mandargs %in% passed)) {
        rlang::abort(paste(
            "missing values for",
            paste(setdiff(mandargs, passed), collapse = ", ")
        ))
    }

    txs <- object[["transcript"]]$transcript_id
    # check input counts features
    if(is.null(rownames(countData))){
        rlang::abort("Input countData do not have feature names")
    } else if(!all(txs %in% rownames(countData))){
        # TODO: consider returning warning instead. But how to handle no exp counts?
        rlang::abort("Some transcript features are missing in countData")
    }

    if(verbose){ .msgheader("Adding expression data")}
    # convert counts to integer
    if(verbose){ .msgsubinfo("Adding transcript counts")}
    countData <- as.matrix(countData)
    counts.names <- rownames(countData)
    counts.samples <- colnames(countData)
    countData <- matrix(as.integer(countData), ncol = length(counts.samples))
    colnames(countData) <- counts.samples
    rownames(countData) <- counts.names

    object@sets$transcript@counts <- countData[txs,]
    object <- .addSampleData(object, sampleData)
    object <- .processCounts(object, verbose)

    return(object)
})


.addSampleData <- function(object, sampleData) {

    object@colData <- data.frame(row.names = rownames(sampleData))
    object@colData$proj.ident <- object@project
    object@active.ident <- "proj.ident"

    # convert strings to factors
    coltypes <- sapply(sampleData, class)
    toconvert <- names(coltypes[coltypes == "character"])
    sampleData[toconvert] <- lapply(sampleData[toconvert],
                                    function(x) factor(x, levels = unique(x)))

    # add sampleData
    object@colData <- cbind(object@colData, sampleData[rownames(object@colData),])

    return(object)
}


.processCounts <- function(object, verbose = FALSE) {

    txcounts <- object@sets$transcript@counts
    # get gene counts
    if(verbose){ .msgsubinfo("Adding gene counts")}
    txdata <- object[["transcript"]]
    genecounts <- rowsum(txcounts[txdata$transcript_id, ], group = txdata$gene_id)
    genecounts <- genecounts[object[["gene"]]$gene_id,]
    object@sets$gene@counts <- genecounts


    # get psi values
    if(verbose){ .msgsubinfo("Adding spliced-event counts")}
    ## normalize tx reads by tx length
    txlengths <- txdata$width
    txcounts.norm <- apply(txcounts, 2, function(x) x/txlengths)

    ## get list of transcripts that are skipped and included
    gtf <- object@transcriptome
    AS.txs.included <- gtf %>% as.data.frame() %>%
        dplyr::filter(type == "AS") %>%
        dplyr::distinct(AS_id, transcript_id)

    ASevents <- granges(object, set = "AS")
    transcripts <- granges(object, set = "transcript")
    transcripts <- transcripts[transcripts$type == "transcript"]

    AS.txs.overlap <- as.data.frame(IRanges::mergeByOverlaps(transcripts, ASevents))
    AS.txs.overlap <- dplyr::distinct(as.data.frame(AS.txs.overlap[c("transcript_id", "AS_id.1")]))

    AS.txs.included.normcounts <- txcounts.norm[AS.txs.included$transcript_id,]
    AS.txs.included.normcounts <- rowsum(AS.txs.included.normcounts, AS.txs.included$AS_id)

    AS.txs.all.normcounts <- txcounts.norm[AS.txs.overlap$transcript_id,]
    AS.txs.all.normcounts <- rowsum(AS.txs.all.normcounts, AS.txs.overlap$AS_id)

    AS.psi <- AS.txs.included.normcounts/AS.txs.all.normcounts
    AS.psi <- signif(AS.psi, digits = 3)
    object@sets$AS@data <- AS.psi[rownames(object@sets$AS@rowData),]
    object@sets$AS@misc$inc <- AS.txs.included.normcounts  #included norm counts
    object@sets$AS@misc$total <- AS.txs.all.normcounts  #included norm counts

    # normalize gene and tx data
    if(verbose){ .msgsubinfo("Normalizing counts")}
    gene.dds <- DESeq2::DESeqDataSetFromMatrix(genecounts, object@colData, ~1)
    gene.dds <- DESeq2::estimateSizeFactors(gene.dds)
    object@sets$gene@data <- DESeq2::counts(gene.dds, normalized = T)

    tx.dds <- DESeq2::DESeqDataSetFromMatrix(txcounts, object@colData, ~1)
    tx.dds <- DESeq2::estimateSizeFactors(tx.dds)
    object@sets$transcript@data <- DESeq2::counts(tx.dds, normalized = T)



    return(object)
}
























