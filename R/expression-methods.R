# TODO: Need to add sample data and refine this
setGeneric("addTxCounts", function(object, countData, sampleData,
                                   verbose = TRUE) standardGeneric("addTxCounts"))
#' Expression function
#'
#' @description
#' A set of functions to incorporate expression data into factR object and
#' tests regulatory potential of alternative splicing events.
#'
#' @details
#' `addTxCounts` converts transcript-level expression counts into gene-level
#' expression counts. If "psiData" is not provided, this function will determine
#' the inclusion level of alternative splicing events (psi) for each sample
#' by calculating the proportion of transcripts containing the splicing event
#' over all transcripts that spans that splicing event. If "psiData" is provided,
#' will check for common exon coordinates and non-overlapping ones detected
#' by factR2 will be set to NA. Row names of "sampleData" dataframe can be unnamed,
#' and this function will search for column variable corresponding to sample names.
#' If unsuccessful or if row names of "sampleData" does not match to column names
#' of "countData", an error will be returned.
#'
#'
#'
#' @param object factR object
#' @param countData Matrix object containing transcript-level expression counts data.
#' @param sampleData Dataframe containing samples information. Dataframe rows
#' can be named to match the sample names in `countData`. If `sampleData` has no
#' row names, function will attempt to pick the column containing sample names
#' and assign it as rownames.
#' @param verbose Boolean value as to whether messages should be printed (Default: TRUE)
#'
#' @rdname factR-exp
#' @name factR-exp
#'
#' @return factRObject with updated counts data and samples metadata.
#'
#' @seealso \code{\link{runfactR}}
#' @export
#'
#' @examples
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
    object <- .addSampleData(object, sampleData, counts.samples)
    object <- .processCounts(object, verbose)

    return(object)
})


.addSampleData <- function(object, sampleData, samples) {

    # check rownames
    object@colData <- data.frame(row.names = samples)
    object@colData$proj.ident <- object@project
    object@active.ident <- "proj.ident"

    # order sampleData according to samples
    if(!all(rownames(sampleData) %in% object@colData)){
        samples.col <- sapply(sampleData, function(x){all(x==samples)})
        if(!any(samples.col)){
            rlang::abort("
            Unable to match sample names to samples metadata.
            Please provide `sampleData` with row.names that
            matches sample names in expression data")
        }
        row.names(sampleData) <- sampleData[[names(samples.col)[samples.col]]]
    }

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
























