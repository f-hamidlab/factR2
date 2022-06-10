setMethod("processCounts", "factR", function(object, verbose = FALSE) {

    # check if count data is present
    if(nrow(object@sets$transcript@counts) == 0){
        rlang::abort("Count data not found in object")
    }

    txcounts <- object@sets$transcript@counts
    # get gene counts
    txdata <- featureData(object, set = "transcript")
    genecounts <- rowsum(txcounts[txdata$transcript_id, ], group = txdata$gene_id)
    genecounts <- genecounts[featureData(object, set = "gene")$gene_id,]
    object@sets$gene@counts <- genecounts


    # get psi values
    ## normalize tx reads by tx length
    txlengths <- txdata$width
    txcounts.norm <- apply(txcounts, 2, function(x) x/txlengths)

    ## get list of transcripts that are skipped and included
    gtf <- object@transcriptome
    AS.txs.included <- gtf %>% as.data.frame() %>%
        dplyr::filter(type == "AS") %>%
        dplyr::mutate(id = paste0(seqnames, ":", start, "-", end, ":", AStype)) %>%
        dplyr::distinct(id, transcript_id)
    AS.txs.skipped <- gtf %>% as.data.frame() %>%
        dplyr::filter(type == "AS") %>%
        dplyr::mutate(id = paste0(seqnames, ":", start, "-", end, ":", AStype)) %>%
        dplyr::distinct(id, gene_id) %>%
        dplyr::left_join(txdata, by = "gene_id") %>%
        dplyr::select(transcript_id, id)
    AS.txs.skipped <- dplyr::setdiff(AS.txs.skipped, AS.txs.included)

    AS.txs.included.normcounts <- txcounts.norm[AS.txs.included$transcript_id,]
    AS.txs.included.normcounts <- rowsum(AS.txs.included.normcounts, AS.txs.included$id)

    AS.txs.skipped.normcounts <- txcounts.norm[AS.txs.skipped$transcript_id,]
    AS.txs.skipped.normcounts <- rowsum(AS.txs.skipped.normcounts, AS.txs.skipped$id)

    AS.txs.sum.normcounts <- AS.txs.included.normcounts + AS.txs.skipped.normcounts
    AS.psi <- AS.txs.included.normcounts/AS.txs.sum.normcounts
    AS.psi <- signif(AS.psi, digits = 3)
    object@sets$AS@data <- AS.psi

    # normalize gene and tx data
    gene.dds <- DESeq2::DESeqDataSetFromMatrix(genecounts, object@colData, object@design)
    gene.dds <- DESeq2::estimateSizeFactors(gene.dds)
    object@sets$gene@data <- DESeq2::counts(gene.dds, normalized = T)

    tx.dds <- DESeq2::DESeqDataSetFromMatrix(txcounts, object@colData, object@design)
    tx.dds <- DESeq2::estimateSizeFactors(tx.dds)
    object@sets$transcript@data <- DESeq2::counts(tx.dds, normalized = T)



    return(object)
})




setMethod("addCountData", "factR", function(object, countData, sampleData, design, verbose = FALSE) {


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
    object <- addDesign(object, design) 

    object <- .processCounts(object) 

    return(object)
})


.addSampleData <- function(object, sampleData) {

    # convert strings to factors 
    data.names <- rownames(sampleData)
    sampleData <- as.data.frame(unclass(testSamples),stringsAsFactors=TRUE)
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
        samples(object) <- cbind(object@colData, sampleData[current_samples,])
    }
    return(object)
}


setMethod("design", "factR", function(object, countData, sampleData, design, verbose = FALSE) {})

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




























