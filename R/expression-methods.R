setGeneric("addTxCounts", function(object, countData, sampleData=NULL, psi=NULL,
                                   verbose = TRUE) standardGeneric("addTxCounts"))
#' Expression/samples-related functions
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
#' @param object factR object
#' @param countData Can be one of the following:
#' \itemize{
#'  \item{}{Path to local count matrix file in .tsv or .csv format}
#'  \item{}{Matrix object containing transcript-level expression counts data.}
#' }
#' @param sampleData (Optional) Dataframe containing samples information. Dataframe rows
#' can be named to match the sample names in `countData`. If `sampleData` has no
#' row names, function will attempt to pick the column containing sample names
#' and assign it as rownames.
#' @param psi (Optional) Exon-level splicing inclusion data. Can be one of the following:
#' \itemize{
#'  \item{}{Path to local file in .tsv or .csv format}
#'  \item{}{Matrix object}
#' }
#' Rownames of matrix should be in chr:start-end format
#' @param verbose Boolean value as to whether messages should be printed (Default: TRUE)
#'
#' @rdname factR-exp
#' @name factR-exp
#'
#' @return factRObject with updated counts data and samples metadata.
#'
#' @seealso \code{\link{factRObject-class}} \code{\link{factR-exp-meta}}
#' @export
#'
#' @examples
#' ## get path to sample GTF and expression data
#' gtf <- system.file("extdata", "pb_custom.gtf.gz", package = "factR2")
#' counts <- system.file("extdata", "pb_expression.tsv.gz", package = "factR2")
#'
#' ## create factRObject with expression counts
#' factRObject <- createfactRObject(gtf, "vM33", countData = counts)
#'
#' ### This can also be done post-creation of a factR object
#' \dontrun{
#' factRObject <- createfactRObject(gtf, "vM33")
#' factRObject <- addTxCounts(factRObject, counts)}
#'
#' ## access counts data
#' counts(factRObject)  # returns normalised expression data from current Set
#' counts(factRObject, "Ptbp1")  # same as above, but only for specific genes
#' counts(factRObject, "Ptbp1", set = "transcript")  # get transcript-level expression
#' counts(factRObject, "Ptbp1", set = "gene") # get gene-level expression
#' counts(factRObject, "Ptbp1", set = "gene", slot = "counts") # get gene-level count
#'
#' ## access samples metadata
#' samples(factRObject)
#' ident(factRObject)   # prints out the current identity
#'
#' ## run correlation between gene expression and exon inclusion
#' factRObject <- testGeneCorr(factRObject)
#' ase(factRObject)
#'
setMethod("addTxCounts", "factR", function(object, countData,
                                           sampleData=NULL,
                                           psi=NULL,
                                           verbose = TRUE) {

    # catch missing args
    mandargs <- c("countData")
    passed <- names(as.list(match.call())[-1])
    if (any(!mandargs %in% passed)) {
        rlang::abort(paste(
            "missing values for",
            paste(setdiff(mandargs, passed), collapse = ", ")
        ))
    }

    object <- .expimport(object, countData, verbose)
    object <- .addSampleData(object, sampleData, verbose)
    object <- .processCounts(object, verbose)
    if(!is.null(psi)){
        object <- .overwritepsi(object,psi, verbose)
    }

    return(object)
})

setGeneric("testGeneCorr", function(
        object,
        vst = TRUE,
        min_n = 3,
        ...) standardGeneric("testGeneCorr"))
#' Assess regulatory potential of AS events
#'
#' @description Correlates exon inclusion levels with gene expression levels
#'
#'
#' @param object factR object
#' @param vst whether to apply variance stabilization on splicing and expression levels
#' @param min_n minimum number of non-NA samples required for correlation testing
#' @param ... additional arguments parsed to cor.test function
#'
#' @return factRObject with updated ASE metadata
#' @export
#' @seealso \code{\link{cor.test}} \code{\link{factR-exp}}
#'
#' @rdname factR-exp
setMethod("testGeneCorr", "factR", function(
        object,
        vst = TRUE,
        min_n = 3,
        ...) {

    return(.ASgenecorr(object, vst, ...))
})


.ASgenecorr <- function(object, vst = TRUE, min_n=3, ...){

    # TODO: test if all genes/events are in object
    psi <- object@sets$AS@data
    n.psi.NA <- rowSums(!is.na(psi))
    passed.ASevents <- names(n.psi.NA[n.psi.NA >= min_n])

    normexp <- object@sets$gene@data

    # get AS-gene match
    AS2gene <- object[["AS"]] %>%
        dplyr::select(AS_id, gene_id, ASNMDtype) %>%
        dplyr::filter(AS_id %in% passed.ASevents)
    psi <- psi[AS2gene$AS_id,]

    if(vst){
        .msgheader("Creating samples metadata")
        psi <- suppressWarnings(.asinTransform(psi))
        ASNMD.stim <- AS2gene$ASNMDtype == "Stimulating"
        ASNMD.stim <- tidyr::replace_na(ASNMD.stim, FALSE)
        psi[ASNMD.stim,] <- 1-psi[ASNMD.stim,]
        normexp <- DESeq2::varianceStabilizingTransformation( object@sets$gene@counts)
    }

    # get samples
    samples <- rownames(object@colData)

    # run correlation
    .msgheader(stringr::str_glue("Running correlation on {nrow(AS2gene) events}"))
    out <- do.call(rbind, pbapply::pbapply(AS2gene, 1, function(dat){
        AS <- dat[1]
        gene <- dat[2]

        psi.dat <- psi[AS,samples]
        exp.dat <- normexp[gene, samples]
        exp.dat <- ifelse(is.na(psi.dat), NA, exp.dat) # match NA

        test <- suppressWarnings(cor.test(psi.dat, exp.dat, ...))
        data.frame(estimate = test$estimate,
                   pvalue = test$p.value)

    }))

    # update ASE df
    object@sets$AS@rowData$gene.cor.estimate <- NA
    object@sets$AS@rowData$gene.cor.pval <- NA
    object@sets$AS@rowData[AS2gene$AS_id, "gene.cor.estimate"] <- out$estimate
    object@sets$AS@rowData[AS2gene$AS_id, "gene.cor.pval"] <- out$pvalue

    object
}




.addSampleData <- function(object, sampleData=NULL, verbose) {

    if(verbose){ .msgheader("Creating samples metadata")}
    # check rownames
    samples <- colnames(object@sets$transcript@counts)
    object@colData <- data.frame(row.names = samples)
    object@colData$proj.ident <- object@project
    object@active.ident <- "proj.ident"

    if(!is.null(sampleData)){
        if(verbose){ .msgsubinfo("Adding provided samples information")}
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
    }
    return(object)
}


.processCounts <- function(object, verbose = FALSE) {

    if(verbose){ .msgheader("Processing expression data")}

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


.expimport <- function(object, countData, verbose){

    if(verbose){ .msgheader("Adding expression data")}


    # check if a path is given
    countDat <- NULL
    if("matrix" %in% class(countData)){
        if(verbose){ .msgsubinfo("Using matrix object")}
        countDat <- as.matrix(countData)
    } else if("character" %in% class(countData)) {
        if(file.exists(countData)){
            if(verbose){ .msgsubinfo("Importing from local file")}
            # check for tsv or csv
            ext <- tools::file_ext(stringr::str_remove(countData, ".gz"))
            if(ext=="tsv"){
                countDat <- as.matrix(utils::read.delim(countData))
            } else if(ext == "csv"){
                countDat <- as.matrix(read.csv(countData))
            }
        }
    }
    if(is.null(countDat)){
        rlang::abort("Input to `countData` not recognized or file does not exist")
    }


    txs <- object[["transcript"]]$transcript_id
    # check input counts features
    if(is.null(rownames(countDat))){
        rlang::abort("Input countData do not have feature names")
    }
    if(!all(txs %in% rownames(countDat))){
        rlang::warn("Some transcript features are missing in countData. Setting expression to '0'")

        # get txs missing from countData
        missing.txs <- txs[!txs %in% rownames(countDat)]
        missing.countData <- matrix(0,
                                    nrow = length(missing.txs),
                                    ncol = ncol(countDat))
        rownames(missing.countData) <- missing.txs
        colnames(missing.countData) <- colnames(countDat)
        countDat <- rbind(countDat, missing.countData)
        countDat <- countDat[txs,]

    }


    # convert counts to integer
    countDat <- as.matrix(countDat)
    counts.names <- rownames(countDat)
    counts.samples <- colnames(countDat)
    countDat <- matrix(as.integer(countDat), ncol = length(counts.samples))
    colnames(countDat) <- counts.samples
    rownames(countDat) <- counts.names

    object@sets$transcript@counts <- countDat[txs,]
    object
}



.overwritepsi  <- function(object, countData, verbose){

    if(verbose){ .msgheader("Using user-defined splicing data")}

    # check if a path is given
    countDat <- NULL
    if("matrix" %in% class(countData)){
        if(verbose){ .msgsubinfo("Using matrix object")}
        countDat <- as.matrix(countData)
    } else if("character" %in% class(countData)) {
        if(file.exists(countData)){
            if(verbose){ .msgsubinfo("Importing from local file")}
            # check for tsv or csv
            ext <- tools::file_ext(stringr::str_remove(countData, ".gz"))
            if(ext=="tsv"){
                countDat <- as.matrix(utils::read.delim(countData))
            } else if(ext == "csv"){
                countDat <- as.matrix(read.csv(countData))
            }
        }
    }
    if(is.null(countDat)){
        rlang::abort("Input to `psi` not recognized or file does not exist")
    }

    as <- object[["AS"]]$coord
    # check input counts features
    if(is.null(rownames(countDat))){
        rlang::abort("Input `psi` do not have feature names")
    }
    if(!all(as %in% rownames(countDat))){
        rlang::warn("Some alternative exon features are missing in psi. Setting PSI to '0'")

        # get txs missing from countData
        missing.txs <- as[!as %in% rownames(countDat)]
        missing.countData <- matrix(0,
                                    nrow = length(missing.txs),
                                    ncol = ncol(countDat))
        rownames(missing.countData) <- missing.txs
        colnames(missing.countData) <- colnames(countDat)
        countDat <- rbind(countDat, missing.countData)
        countDat <- countDat[as,]

    }


    # convert counts to integer
    countDat <- as.matrix(countDat)
    rownames(countDat) <- object[["AS"]]$AS_id

    object@sets$AS@data <- countDat
    object
}

















