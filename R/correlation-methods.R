setGeneric("testGeneCorr", function(object, ...) standardGeneric("testGeneCorr"))
setMethod("testGeneCorr", "factR", function(
        object,
        vst = TRUE,
        ...) {

    return(.ASgenecorr(object, vst, ...))
})


.ASgenecorr <- function(object, vst = TRUE, ...){

    # get AS-gene match
    AS2gene <- object[["AS"]] %>%
        dplyr::select(AS_id, gene_id)

    psi <- object@sets$AS@data
    normexp <- object@sets$gene@data
    if(vst){
        psi <- .asinTransform(psi)
        normexp <- DESeq2::varianceStabilizingTransformation( object@sets$gene@counts)
    }

    # get samples
    samples <- rownames(object@colData)

    # run correlation
    out <- do.call(rbind, apply(AS2gene, 1, function(dat){
        AS <- dat[1]
        gene <- dat[2]

        test <- cor.test(psi[AS,samples], normexp[gene, samples], ...)
        data.frame(estimate = test$estimate,
                   pvalue = test$p.value)

    }))

    # update ASE df
    object@sets$AS@rowData$gene.cor.estimate <- out$estimate
    object@sets$AS@rowData$gene.cor.pval <- out$pvalue

    object
}
