## TODO: make 1-asin a contextual arg

#' Assess regulatory potential of AS events
#'
#' @description Correlates inclusion levels of exons with gene expression levels
#'
#'
#' @param object factR object
#' @param vst whether to apply variance stabilization on splicing and expression levels
#' @param min_n minimum number of non-NA samples required for correlation testing
#' @param ... additional arguments parsed to cor.test function
#'
#' @return factRObject with updated ASE metadata
#' @export
#' @seealso \code{\link{cor.test}}
#'
#' @rdname testGeneCorr
#' @examples
#' data(factRsample)
#' factRsample <- runfactR(factRsample)
setGeneric("testGeneCorr", function(
        object,
        vst = TRUE,
        min_n = 3,
        ...) standardGeneric("testGeneCorr"))
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
        psi <- suppressWarnings(.asinTransform(psi))
        ASNMD.stim <- AS2gene$ASNMDtype == "Stimulating"
        ASNMD.stim <- tidyr::replace_na(ASNMD.stim, FALSE)
        psi[ASNMD.stim,] <- 1-psi[ASNMD.stim,]
        normexp <- DESeq2::varianceStabilizingTransformation( object@sets$gene@counts)
    }

    # get samples
    samples <- rownames(object@colData)

    # run correlation
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
