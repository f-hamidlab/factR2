setGeneric("getAAsequence", function(object, verbose = FALSE) standardGeneric("getAAsequence"))
setMethod("getAAsequence", "factR", function(object, verbose = FALSE) {
    gtf <- object@transcriptome
    if(! "CDS" %in% gtf$type){
        rlang::abort("No CDSs found. Please run buildCDS() first")
    }

    genetxs <- object[["transcript"]]
    txs <- genetxs[genetxs$cds == "yes",]$transcript_id

    gtf <- gtf[gtf$transcript_id %in% txs]
    if("DNAStringSet" %in% is(object@reference$genome)){
        chrnames <- names(object@reference$genome)
    } else {
        chrnames <- GenomeInfoDb::seqnames(object@reference$genome)
    }

    gtf <- gtf[as.character(GenomeInfoDb::seqnames(gtf)) %in% chrnames]
    cds <- S4Vectors::split(gtf[gtf$type == "CDS"], ~transcript_id)



    out <- .getSequence(cds, object@reference$genome, verbose)
    outseq <- Biostrings::AAStringSet(out$x)
    names(outseq) <- out$id
    object@domains$sequence <- outseq

    return(object)
})


setGeneric("predictDomain", function(object, ...,
                                     database = "superfamily",
                                     ncores = 4) standardGeneric("predictDomain"))
setMethod("predictDomain", "factR", function(object,
                                              ...,
                                              database = "superfamily",
                                              ncores = 4) {

    # check if CDS have been built
    gtf <- object@transcriptome
    if(! "CDS" %in% gtf$type){
        rlang::abort("No CDSs found. Please run buildCDS() first")
    }

    # get transcripts to test
    genetxs <- features(object, ..., set = "transcript")
    txs <- genetxs$transcript_id


    # see if txs to test have already been tested
    tested <- object@domains$tested
    txs.to.test <- txs[!txs %in% tested]
    object@domains$tested <- c(tested, txs.to.test)

    if(length(txs.to.test) > 0){
        cdss <- object@domains$sequence
        cdss.to.test <- cdss[which(names(cdss) %in% txs.to.test)]

        if(length(cdss.to.test) == 0){
            rlang::warn("None of the transcripts to test contain CDSs")
            return(object)
        } else if(length(cdss.to.test) < length(txs.to.test)){
            rlang::warn(sprintf("Skipped %s non-coding transcripts",
                                length(txs.to.test)-length(cdss.to.test)))
        }


        if(length(cdss.to.test) > 1000){
            continue <- ""
            while(!continue %in% c("y","n")){
                continue <- readline("Predicting domains for >1000 proteins. Continue? [y/n]  ")
            }
            if(continue == "n"){
                rlang::abort("Function aborted")
            }
        }
        # run domain prediction
        domains.out <- .runDomainSearch(cdss.to.test, database, ncores)
        domains <- slot(object, "domains")$data
        slot(object, "domains")$data <- rbind(domains, domains.out)
    }
    object
})




.getSequence <- function(cds, fasta, verbose = FALSE) {
    x <- y <- instop <- NULL

    if(verbose){
        .msgheader("Translating amino acid sequences")
    }

    cdsSeq <- GenomicFeatures::extractTranscriptSeqs(fasta, cds)
    aaSeq <- suppressWarnings(
        Biostrings::translate(cdsSeq, if.fuzzy.codon = "solve")) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("id") %>%
        dplyr::rowwise() %>%
        dplyr::mutate(y = strsplit(x, split = "")) %>%
        dplyr::mutate(noATG = ifelse(y[[1]] != "M", TRUE, FALSE)) %>%
        dplyr::mutate(instop = ifelse("*" %in% y, TRUE, FALSE)) %>%
        dplyr::ungroup()

    # check for ATG and internal stop_codon, truncate proteins with internal
    # stop codon
    ## and remove entries without proteins after truncation
    if (TRUE %in% aaSeq$noATG) {
        if(verbose){ .msgsubwarn(
            sprintf("%s CDSs do not begin with ATG", sum(aaSeq$noATG)))}
    }
    if (TRUE %in% aaSeq$instop) {
        aaSeq <- suppressWarnings(aaSeq %>%
                                      dplyr::rowwise() %>%
                                      dplyr::mutate(x = ifelse(instop == TRUE,
                                                               paste(y[seq_len(which(y == "*") - 1)], collapse = ""),
                                                               x
                                      )) %>%
                                      dplyr::mutate(y = strsplit(x, split = "")) %>%
                                      dplyr::ungroup())
        if(verbose){ .msgsubwarn(
            sprintf(paste0("%s CDSs contain internal stop codon. ",
                           "Truncating CDS sequence to retain ORF"),
                    sum(aaSeq$instop)))}

        if ("" %in% aaSeq$x) {
            if(verbose){ .msgsubwarn(
                sprintf(paste0(
                    "After truncation, %s cds have no ",
                    "coding sequences. These CDSs were not analyzed"),
                    sum(aaSeq$x == "")))}

            aaSeq <- aaSeq[aaSeq$x != "", ]
        }
    }

    return(aaSeq)
}



.runDomainSearch <- function(aaSeq, db = "superfamily", numcores = 1) {
    type <- transcript_id <- description <- begin <- id <- NULL

    # prepare URL
    url <- paste("https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan")
    url.opts <- list(
        httpheader = "Expect:", httpheader = "Accept:text/xml",
        verbose = FALSE, followlocation = TRUE
    )

    # run search for each protein sequence
    output <- BiocParallel::bplapply(names(aaSeq), function(y) {
        #print(aaSeq[y,]$id)

        # account for return errors
        report <- tryCatch(
            .getdomains(url, url.opts, as.character(aaSeq[[y]]),
                        y, db),
            error = function(e) NULL
        )

        if (is.null(report)) {
            return(tibble::tibble(type = "CHAIN",
                                  description = y,
                                  eval = NA,
                                  begin = 1,
                                  end = length(aaSeq[[y]]),
                                  transcript_id = y))
        } else {
            return(dplyr::bind_rows(
                report,
                tibble::tibble(type = "CHAIN",
                               description = y,
                               eval = NA,
                               begin = 1,
                               end = length(aaSeq[[y]]),
                               transcript_id = y)
            ))
        }
    }, BPPARAM = BiocParallel::MulticoreParam(tasks = numcores)) %>%
        dplyr::bind_rows()

    return(output)
}


.getdomains <- function(url, url.opts, seq, id, db = "superfamily") {
    type <- famdesc <- fameval <- begin <- NULL

    hmm <- ""
    while(!grepl("data name=\"results\"", hmm)){
        hmm <- RCurl::postForm(url,
                               hmmdb = db, seqdb = NULL,
                               seq = seq, style = "POST", .opts = url.opts,
                               .contentEncodeFun = RCurl::curlPercentEncode,
                               .checkParams = TRUE
        )
        if (grepl("status=\"PEND\"", hmm)) {
            Sys.sleep(30)
        }
    }


    xml <- XML::xmlParse(gsub("<act_site.*?/act_site>\n", "", hmm))
    domains <- XML::xpathSApply(xml, "///domains", XML::xpathSApply, "@*")
    if("matrix" %in% class(domains)){
        vector <- as.character(domains)
        names(vector) <- rownames(domains)
        domains <- list(vector)
    }
    purrr::map_dfr(domains, function(x){
        x <- as.list(x)
        return(data.frame(type = "DOMAIN",
                          description = x$alihmmdesc,
                          eval = x$cevalue,
                          begin = as.numeric(x$alisqfrom),
                          end = as.numeric(x$alisqto),
                          transcript_id = id
                          ))
    })
}
