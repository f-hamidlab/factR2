createfactRObject <- function(gtf, reference,
                              use_own_annotation = NULL,
                              use_own_genome = NULL,
                              countData = NULL,
                              sampleData = NULL,
                              design = NULL,
                              psi = NULL,
                              verbose = FALSE){

    # check input arguments
    if(verbose){
        rlang::inform("Checking inputs")
    }
    ## GTF input
    if(missing(gtf)){
        rlang::abort("Missing gtf input")
    }
    ## Annotation input
    reference.list <- .checkannotation(reference, use_own_annotation,
                                       use_own_genome, verbose)

    ## Count data input
    if(!is.null(countData)){
        # check if sampleData and design are given
        if(any(c(is.null(sampleData), is.null(design)))){
            args <- c("sampleData", "design")
            args.null <- c(is.null(sampleData), is.null(design))
            rlang::abort(sprintf("Counts data provided, but %s data missing",
                                 paste(args[args.null], collapse = " and ")))
        }
    }


    # create new factRObject
    if(verbose){
        rlang::inform("Creating factRObject")
    }
    obj <- methods::new("factR")
    obj@version <- factR2version

    # import data
    if(verbose){
        rlang::inform("Adding custom transcriptome")
    }
    obj@transcriptome <- .smartimport(gtf, ".gtf", "Custom transcriptome",verbose)

    if(verbose){
        rlang::inform("Adding annotation")
    }
    obj@reference$ranges <- .smartimport(reference.list[[1]], ".gtf", "Annotation",
                                         verbose)
    if(verbose){
        rlang::inform("Adding genome sequence")
    }
    obj@reference$genome <- .smartimport(reference.list[[2]], ".fa", "Genome sequence",
                                         verbose)


    # run object prepration function
    obj <- .prepfactR(obj, verbose)

    # add and prep counts data if given
    if(!is.null(countData)){
        if(verbose){ rlang::inform("Adding expression counts data")}
        obj <- addTxCounts(obj, countData, sampleData, design)
    }

    return(obj)
}

.checkannotation <- function(ref, own_annotation, own_genome, verbose){
    annotation <- ""
    genome <- ""
    if(missing(ref)){
        empty_custom_ref <- c(is.null(own_annotation), is.null(own_genome))
        if(all(empty_custom_ref)){
            rlang::abort("Missing value for reference")
        } else if(any(empty_custom_ref)){
            rlang::abort(sprintf("Using own reference, but value for '%s' not provided",
                                 c("use_own_annotation", "use_own_genome")[empty_custom_ref]))
        } else {
            annotation <- own_annotation
            genome <- own_genome
        }
    }
    else {
        data(genomes)
        # select genome based on given reference
        if(ref %in% genomes$ID){
            selected_genome <- genomes[genomes$ID == ref,]
        } else if(any(stringr::str_detect(genomes$synonyms, ref))){
            selected_genome_index <- which( stringr::str_detect(genomes$synonyms,
                                                                ref))
            selected_genome <- genomes[selected_genome_index,]
        } else {
            rlang::abort(sprintf("Input reference '%s' not recognized",
                                 ref))
        }

        annotation <- selected_genome$annotation
        ## test if BSgenome object is installed
        if(selected_genome$genome.pri %in% BSgenome::installed.genomes()){
            require(selected_genome$genome.pri, character.only = T,
                    quietly = !verbose)
            genome <- get(selected_genome$genome.pri)
        } else {
            genome <- selected_genome$genome.sec
        }
    }
    return(list(annotation, genome))

}



.smartimport <- function(input, format, type, verbose=FALSE){
    options(timeout = 10000)

    # import as an object
    if(is.object(input)){
        if(format == ".gtf"){

            #check if object is a GTF
            if(!is_gtf(input)){
                rlang::abort(sprintf("%s is not a GRanges GTF object", type))
            }
        }
        if(verbose){
            datatype <- as.character(class(input))
            rlang::inform(sprintf("## Using %s object", datatype))
        }
        return(input)
    }



    ## try local file
    else if (file.exists(input) & grepl(format, tolower(input))){
        if(verbose){
            rlang::inform("## Importing from local directory")
        }
        if(format == ".gtf"){
            return(factR::importGTF(input))
        } else if(format == ".fa"){
            return(factR::importFASTA(input))
        }

    }
    ## try AH
    else if(grepl("^AH", input)){
        if(verbose){
            rlang::inform("## Importing from AnnotationHub")

            ah <- AnnotationHub::AnnotationHub()
            return(ah[[input]])
        } else {
            ah <- suppressMessages(AnnotationHub::AnnotationHub())
            return(suppressMessages(ah[[input]]))
        }

    }
    ## try FTP
    else if(grepl("^http",input)){
        if(verbose){
            rlang::inform("## Importing from URL")
        }
        if(format == ".gtf"){
            tmpfile <- tempfile()
            download.file(input, tmpfile, quiet = !verbose)
            return(factR::importGTF(tmpfile))
        } else if(format == ".fa"){
            tmpfile <- tempfile()
            download.file(input, tmpfile, quiet = !verbose)
            return(factR::importFASTA(tmpfile))
        }

    } else {
        rlang::abort(sprintf("Input not recognized or file does not exists"))
    }
}



.prepfactR <-  function(object, verbose = FALSE) {
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
    if("ref_gene_id" %in% colnames(as.data.frame(object@transcriptome))){
        if(verbose){
            object@transcriptome <- factR::matchGeneInfo(object@transcriptome,
                                                         object@reference$ranges,
                                                         primary_gene_id = "gene_id",
                                                         secondary_gene_id = "ref_gene_id")
        } else {
            object@transcriptome <- suppressMessages(
                factR::matchGeneInfo(object@transcriptome,
                                     object@reference$ranges,
                                     primary_gene_id = "gene_id",
                                     secondary_gene_id = "ref_gene_id"))
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





    # custom.df <- as.data.frame(object@transcriptome)
    # potential_id_vars <- apply(custom.df, 2, function(x) any(grepl("ENS",x)))
    # potential_id <- names(potential_id_vars)[potential_id_vars]
    # potential_id <- potential_id[-which( potential_id %in% "transcript_id")]


    # if(length(potential_id > 1)){
    #     if(verbose){
    #         object@transcriptome <- factR::matchGeneInfo(object@transcriptome,
    #                                                      object@reference$ranges,
    #                                                      primary_gene_id = "gene_id",
    #                                                      secondary_gene_id = potential_id)
    #     } else {
    #         object@transcriptome <- suppressMessages(
    #             factR::matchGeneInfo(object@transcriptome,
    #                                  object@reference$ranges,
    #                                  primary_gene_id = "gene_id",
    #                                  secondary_gene_id = potential_id))
    #     }
    #
    # } else {
    #     if(verbose){
    #         object@transcriptome <- factR::matchGeneInfo(object@transcriptome,
    #                                                      object@reference$ranges,
    #                                                      primary_gene_id = "gene_id")
    #     } else {
    #         object@transcriptome <- suppressMessages(
    #             factR::matchGeneInfo(object@transcriptome,
    #                                  object@reference$ranges,
    #                                  primary_gene_id = "gene_id"))
    #     }
    # }

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
        rlang::inform("Creating factRset objects")
    }
    if(verbose){
        rlang::inform("## Adding gene information")
    }
    object@active.set <- "transcript"
    object@sets$gene <- methods::new("factRset")
    object@sets$gene@rowData <- as.data.frame(object@transcriptome) %>%
        dplyr::filter(type %in% "gene") %>%
        dplyr::select(gene_id, gene_name, strand, width, match_level) %>%
        dplyr::distinct()
    rownames(object@sets$gene@rowData) <- object@sets$gene@rowData$gene_id
    object@sets$gene@counts <- as.matrix(data.frame(row.names =  rownames(object[["gene"]])))
    object@sets$gene@data <- as.matrix(data.frame(row.names =  rownames(object[["gene"]])))

    if(verbose){
        rlang::inform("## Adding transcript information")
    }
    object@sets$transcript <- methods::new("factRset")
    object@sets$transcript@rowData <- as.data.frame(object@transcriptome) %>%
        dplyr::filter(type %in% "exon") %>%
        dplyr::group_by(gene_id, gene_name, transcript_id, strand) %>%
        dplyr::summarise(width = sum(width)) %>%
        dplyr::select(transcript_id, gene_id, gene_name, strand, width) %>%
        as.data.frame()
    rownames(object@sets$transcript@rowData) <- object@sets$transcript@rowData$transcript_id
    object@sets$transcript@counts <- as.matrix(data.frame(row.names =  rownames(object[["transcript"]])))
    object@sets$transcript@data <- as.matrix(data.frame(row.names =  rownames(object[["transcript"]])))

    if(verbose){
        rlang::inform("## Adding alternative splicing information")
    }
    object@sets$AS <- methods::new("factRset")
    object <- .findAS(object)
    object@sets$AS@rowData <- as.data.frame(object@transcriptome) %>%
        dplyr::filter(type %in% "AS") %>%
        dplyr::mutate(coord = paste0(seqnames, ":", start, "-", end)) %>%
        dplyr::select(AS_id, gene_id, gene_name, coord, AStype, strand, width) %>%
        dplyr::distinct() %>%
        dplyr::mutate(AStype = factor(AStype, levels = c("CE", "AD","AA","AF","AL","RI")))
    rownames(object@sets$AS@rowData) <- object@sets$AS@rowData$AS_id
    object@sets$AS@counts <- as.matrix(data.frame(row.names =  rownames(object[["AS"]])))
    object@sets$AS@data <- as.matrix(data.frame(row.names =  rownames(object[["AS"]])))

    # annotate new transcripts
    if(verbose){
        rlang::inform("Annotating novel transcripts")
    }
    gtf <- granges(object, set = "transcript")
    newtxs <- suppressMessages(factR::subsetNewTranscripts(gtf,
                                                           object@reference$ranges,
                                                           refine.by = "intron"))
    object@sets$transcript@rowData$novel <- ifelse(object[[]]$transcript_id %in% newtxs$transcript_id,
                                                   "yes", "no")
    object@sets$transcript@rowData$cds <- "no"
    object@sets$transcript@rowData$nmd <- "no"
    object
}
















