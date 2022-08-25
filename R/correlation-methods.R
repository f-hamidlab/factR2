## TODO: do up documentation

setGeneric("testGeneCorr", function(object, ...) standardGeneric("testGeneCorr"))
setMethod("testGeneCorr", "factR", function(
        object,
        vst = TRUE,
        min_n = 3,
        ...) {

    return(.ASgenecorr(object, vst, ...))
})


.ASgenecorr <- function(object, vst = TRUE, ...){

    psi <- object@sets$AS@data
    n.psi.NA <- rowSums(!is.na(psi))
    passed.ASevents <- names(n.psi.NA[n.psi.NA >= min_n])

    normexp <- object@sets$gene@data

    # get AS-gene match
    AS2gene <- object[["AS"]] %>%
        dplyr::select(AS_id, gene_id) %>%
        dplyr::filter(AS_id %in% passed.ASevents)

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
