setMethod("runfactR", "factR", function(object, verbose = FALSE) {
    object <- prepTranscriptome(object, verbose)
    object <- buildCDS(object, verbose)
    object <- predictNMD(object, verbose)
    object <- getAAsequence(object, verbose)
    object <- testASNMDevents(object)
    object
})


setMethod("prepTranscriptome", "factR", function(object, verbose = FALSE) {
    # match chromosomes and gene ID
    if(verbose){
        rlang::inform("Matching chromosome names")
    }
    object@transcriptome <- factR::matchChromosomes(object@transcriptome,
                                                 object@reference$genome)
    object@reference$ranges <- suppressWarnings(factR::matchChromosomes(object@reference$ranges,
                                                                     object@reference$ranges))
    ## try find variables that contain gene ids
    if(verbose){
        rlang::inform("Matching gene names")
    }
    custom.df <- as.data.frame(object@transcriptome)
    potential_id_vars <- apply(custom.df, 2, function(x) any(grepl("ENS",x)))
    potential_id <- names(potential_id_vars)[potential_id_vars]
    potential_id <- potential_id[-which("transcript_id" %in% potential_id)]
    if(length(potential_id > 1)){
        if(verbose){
            object@transcriptome <- factR::matchGeneInfo(object@transcriptome,
                                                      object@reference$ranges,
                                                      primary_gene_id = "gene_id",
                                                      secondary_gene_id = potential_id)
        } else {
            object@transcriptome <- suppressMessages(
                factR::matchGeneInfo(object@transcriptome,
                                     object@reference$ranges,
                                     primary_gene_id = "gene_id",
                                     secondary_gene_id = potential_id))
        }

    } else {
        if(verbose){
            object@transcriptome <- factR::matchGeneInfo(object@transcriptome,
                                                      object@reference$ranges,
                                                      primary_gene_id = "gene_id")
        } else {
            object@transcriptome <- suppressMessages(
                factR::matchGeneInfo(object@transcriptome,
                                     object@reference$ranges,
                                     primary_gene_id = "gene_id"))
        }
    }

    # Add 'gene' type in GTF
    if(!"gene" %in% object@transcriptome$type){
        genes.gtf <- object@transcriptome
        gene.names <- genes.gtf %>%
            as.data.frame() %>%
            dplyr::group_by(gene_id, gene_name) %>%
            dplyr::arrange(dplyr::desc(match_level)) %>%
            dplyr::distinct(gene_id, gene_name, .keep_all = T) %>%
            dplyr::ungroup() %>%
            dplyr::select(gene_id, gene_name, match_level)
        genes.list <- range(
            S4Vectors::split(genes.gtf[genes.gtf$type %in% "transcript"],
                             ~gene_id))
        genes.gtf <- genes.list %>%
            as.data.frame() %>%
            dplyr::group_by(group) %>%
            dplyr::mutate(start = min(start), end = max(end)) %>%
            dplyr::ungroup() %>%
            dplyr::rename(gene_id = group_name) %>%
            dplyr::select(-width, -group) %>%
            dplyr::mutate(source = "factR2", type = "gene") %>%
            dplyr::distinct(gene_id, .keep_all = T) %>%
            dplyr::left_join(gene.names, by = "gene_id") %>%
            GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
        object@transcriptome <- c(genes.gtf, object@transcriptome)
    }


    # create genetxs dataframe
    if(verbose){
        rlang::inform("Creating assays")
    }
    object@active.set <- "transcript"
    object@sets$gene <- methods::new("factRset")
    object@sets$gene@rowData <- as.data.frame(object@transcriptome) %>%
        dplyr::filter(type %in% "gene") %>%
        dplyr::select(gene_id, gene_name, width, match_level) %>%
        dplyr::distinct()
    rownames(object@sets$gene@rowData) <- object@sets$gene@rowData$gene_id

    object@sets$transcript <- methods::new("factRset")
    object@sets$transcript@rowData <- as.data.frame(object@transcriptome) %>%
        dplyr::filter(type %in% "transcript") %>%
        dplyr::select(gene_id, gene_name, transcript_id, width) %>%
        dplyr::distinct()
    rownames(object@sets$transcript@rowData) <- object@sets$transcript@rowData$transcript_id

    object@sets$AS <- methods::new("factRset")
    object <- findAltSplicing(object)
    object@sets$AS@rowData <- as.data.frame(object@transcriptome) %>%
        dplyr::filter(type %in% "AS") %>%
        dplyr::mutate(coord = paste0(seqnames, ":", start, "-", end)) %>%
        dplyr::select(gene_id, gene_name, coord, AStype, width) %>%
        dplyr::distinct()
    rownames(object@sets$AS@rowData) <- paste0(object@sets$AS@rowData$coord,
                                               ":",object@sets$AS@rowData$AStype)

    # annotate new transcripts
    newtxs <- suppressMessages(factR::subsetNewTranscripts(object@transcriptome,
                                                           object@reference$ranges,
                                                           refine.by = "intron"))
    featureData(object)$novel <- ifelse(featureData(object)$transcript_id %in% newtxs$transcript_id,
                                        "yes", "no")
    featureData(object)$cds <- "no"
    featureData(object)$nmd <- "no"
    object
})
















