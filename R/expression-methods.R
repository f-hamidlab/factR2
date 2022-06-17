



setGeneric("addTxCounts", function(object, countData, sampleData,
                                   design, verbose = FALSE) standardGeneric("addTxCounts"))
setMethod("addTxCounts", "factR", function(object, countData, sampleData, design, verbose = FALSE) {

    # catch missing args
    mandargs <- c("countData", "sampleData", "design")
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
        rlang::abort("Some transcript features are missing in countData")
    }

    # convert counts to integer
    counts.names <- rownames(countData)
    counts.samples <- colnames(countData)
    countData <- matrix(as.integer(countData), ncol = length(counts.samples))
    colnames(countData) <- counts.samples
    rownames(countData) <- counts.names

    object@sets$transcript@counts <- countData[txs,]
    object@colData <- data.frame(row.names = colnames(countData))
    object <- .addSampleData(object, sampleData)
    object <- .addDesign(object, design)
    object <- .processCounts(object)

    return(object)
})


.addSampleData <- function(object, sampleData) {

    # convert strings to factors
    data.names <- rownames(sampleData)
    sampleData <- as.data.frame(unclass(sampleData),stringsAsFactors=TRUE)
    rownames(sampleData) <- data.names


    # add sampleData if no data is currently present
    if(nrow(samples(object)) == 0){
        samples(object) <- sampleData
    }
    # else, check for sample consistencies
    else {
        current_samples <- rownames(samples(object))

        if(!identical(rownames(sampleData), as.character(c(1:nrow(sampleData))))){
            if(!all(current_samples %in% rownames(sampleData))){
                rlang::abort("Some samples are missing in sampleData")
            }
        } else {
            if(nrow(sampleData) != length(current_samples)){
                rlang::abort(sprintf("Number of samples in sampleData (%s) is not equal to %s",
                                     nrow(sampleData), length(current_samples)))
            } else {
                poss_sample_vars <- apply(sampleData, 2, function(x) all(current_samples %in% x))
                poss_var <- names(poss_sample_vars)[poss_sample_vars]

                if(length(poss_var) > 0){
                    rownames(sampleData) <- sampleData[[poss_var]]

                } else {
                    rlang::abort("Unable to resolve sample names. Please add missing rownames.")
                }
            }
        }
        # actual appending of data
        object@colData <- cbind(object@colData, sampleData[current_samples,])
    }
    return(object)
}



setMethod("design<-", "factR", function(object, value) .addDesign(object, value))

.addDesign <- function(object, design) {
    design <- formula(design)
    sample.meta <- samples(object)

    # check if design variables are in samples
    if(nrow(sample.meta) == 0){
        rlang::abort("Samples metadata not constructed yet")
    } else {
        vars <- as.character(design)
        vars <- strsplit(vars[2], "[[:punct:]]| ")[[1]]
        vars <- vars[vars != ""]

        if(!all(vars %in% colnames(sample.meta))){
            rlang::abort("Some design variables not in samples metadata")
        } else {
            object@design <- design
        }
    }
    return(object)
}



.processCounts <- function(object, verbose = FALSE) {

    txcounts <- object@sets$transcript@counts
    # get gene counts
    txdata <- object[["transcript"]]
    genecounts <- rowsum(txcounts[txdata$transcript_id, ], group = txdata$gene_id)
    genecounts <- genecounts[object[["gene"]]$gene_id,]
    object@sets$gene@counts <- genecounts


    # get psi values
    ## normalize tx reads by tx length
    txlengths <- txdata$width
    txcounts.norm <- apply(txcounts, 2, function(x) x/txlengths)

    ## get list of transcripts that are skipped and included
    gtf <- object@transcriptome
    AS.txs.included <- gtf %>% as.data.frame() %>%
        dplyr::filter(type == "AS") %>%
        dplyr::distinct(ASid, transcript_id)

    ASevents <- granges(object, set = "AS")
    transcripts <- granges(object, set = "transcript")
    transcripts <- transcripts[transcripts$type == "transcript"]

    AS.txs.overlap <- as.data.frame(IRanges::mergeByOverlaps(transcripts, ASevents))
    AS.txs.overlap <- dplyr::distinct(as.data.frame(AS.txs.overlap[c("transcript_id", "ASid.1")]))

    AS.txs.included.normcounts <- txcounts.norm[AS.txs.included$transcript_id,]
    AS.txs.included.normcounts <- rowsum(AS.txs.included.normcounts, AS.txs.included$ASid)

    AS.txs.all.normcounts <- txcounts.norm[AS.txs.overlap$transcript_id,]
    AS.txs.all.normcounts <- rowsum(AS.txs.all.normcounts, AS.txs.overlap$ASid)

    AS.psi <- AS.txs.included.normcounts/AS.txs.all.normcounts
    AS.psi <- signif(AS.psi, digits = 3)
    object@sets$AS@data <- AS.psi[rownames(object@sets$AS@rowData),]

    # normalize gene and tx data
    gene.dds <- DESeq2::DESeqDataSetFromMatrix(genecounts, object@colData, object@design)
    gene.dds <- DESeq2::estimateSizeFactors(gene.dds)
    object@sets$gene@data <- DESeq2::counts(gene.dds, normalized = T)

    tx.dds <- DESeq2::DESeqDataSetFromMatrix(txcounts, object@colData, object@design)
    tx.dds <- DESeq2::estimateSizeFactors(tx.dds)
    object@sets$transcript@data <- DESeq2::counts(tx.dds, normalized = T)



    return(object)
}
























