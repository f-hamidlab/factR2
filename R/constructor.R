#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Construct a factRObject
#'
#' @param gtf path to custom GTF file.
#' @param reference a character value of supported reference ID. Run
#' \code{\link{ListSupportedGenomes}} for a list of supported reference. Input
#' can also be a reference species name ("Mouse" or "Mus musculus or "Mmusculus").
#' To use your own annotation and genome, leave this parameter blank and provide
#' inputs to `use_own_annotation` and `use_own_genome`
#' @param use_own_annotation can be a (1) GenomicRanges object containing
#' reference transcript annotation, (2) path to local GTF file, (3) AnnotationHub
#' ID or (4) URL to an online GTF file.
#' @param use_own_genome can be a (1) BSgenome or DNAStringSet object containing
#' genome sequence, (2) path to local genome fasta file, (3) AnnotationHub
#' ID or (4) URL to an online genome fasta file.
#' @param verbose if TRUE, show progress messages
#'
#' @return factR object
#' @export
#'
#' @examples
#' # Create factRObject using sample GTF and supported reference
#' gtf <- system.file("extdata", "sc_merged_sample.gtf.gz", package = "factR")
#' obj <- createfactRObject(gtf, "vM25")
#'
#' # Create factRObject using custom reference
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' obj <- createfactRObject(gtf, use_own_annotation = "AH49547", use_own_genome = Mmusculus)
createfactRObject <- function(gtf, reference,
                              use_own_annotation = NULL,
                              use_own_genome = NULL,
                              counts,
                              psi,
                              assays = "transcript",
                              verbose = FALSE){

    # check input GTF arguments
    if(missing(gtf)){
        rlang::abort("Missing gtf input")
    }

    # check for reference arguments
    annotation <- ""
    genome <- ""
    if(missing(reference)){
        empty_custom_ref <- c(is.null(use_own_annotation), is.null(use_own_genome))
        if(all(empty_custom_ref)){
            rlang::abort("Missing value for reference")
        } else if(any(empty_custom_ref)){
            rlang::abort(sprintf("Using own reference, but value for '%s' not provided",
                c("use_own_annotation", "use_own_genome")[empty_custom_ref]))
        } else {
            annotation <- use_own_annotation
            genome <- use_own_genome
        }
    }
    else {
        data(genomes)
        # select genome based on given reference
        if(reference %in% genomes$ID){
            selected_genome <- genomes[genomes$ID == reference,]
        } else if(any(stringr::str_detect(genomes$synonyms, reference))){
            selected_genome_index <- which( stringr::str_detect(genomes$synonyms,
                                                                reference))
            selected_genome <- genomes[selected_genome_index,]
        } else {
            rlang::abort(sprintf("Input reference '%s' not recognized",
                                 reference))
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

    # check existence of GTF input and try to import
    if(is.object(gtf)){
        if(!is_gtf(gtf)){
            rlang::abort("Input is not a GRanges GTF object")
        } else if(verbose){
            rlang::inform("Using custom GRanges GTF object")
        }
    }
    else if (!file.exists(gtf)){
        rlang::abort(sprintf("File '%s' does not exist", gtf))
    } else if (!grepl(".gtf", tolower(gtf))){
        rlang::abort(sprintf("File '%s' is not a GTF", gtf))
    } else {
        if(verbose){
            rlang::inform("Importing custom GTF")
        }
        gtf <- factR::importGTF(gtf)
    }
    # create new factRObject
    obj <- methods::new("factR")
    obj@version <- factR2version
    obj@transcriptome <- gtf  # import GTF
    
    
    custom.gtf <- .smartimport(annotation, ".gtf", verbose)
    if(is_gtf(custom.gtf)){
        obj@reference$ranges <- custom.gtf
    } else {
        rlang::abort("Annotation object is not a GTF")
    }
    
    if(verbose){
        rlang::inform("Retrieving reference genome")
    }
    obj@reference$genome <- .smartimport(genome, ".fa",
                                         verbose)
    
    obj <- .prepTranscriptome(obj, verbose)

    return(obj)
}

.smartimport <- function(input, format, verbose=FALSE){
    options(timeout = 10000)

    # import as an object
    if(is.object(input)){
        return(input)
    }
    ## try local file
    else if (file.exists(input) & grepl(format, tolower(input))){
        if(format == ".gtf"){
            return(factR::importGTF(input))
        } else if(format == ".fa"){
            return(factR::importFASTA(input))
        }

    }
    ## try AH
    else if(grepl("^AH", input)){
        if(verbose){
            ah <- AnnotationHub::AnnotationHub()
            return(ah[[input]])
        } else {
            ah <- suppressMessages(AnnotationHub::AnnotationHub())
            return(suppressMessages(ah[[input]]))
        }

    }
    ## try FTP
    else if(grepl("^http",input)){
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
        rlang::abort(sprintf("Input annotation '%s' not recognized",
                             input))
    }
}



.prepTranscriptome <-  function(object, verbose = FALSE) {
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
        rlang::inform("Creating sets")
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
        dplyr::filter(type %in% "exon") %>%
        dplyr::group_by(gene_id, gene_name, transcript_id) %>%
        dplyr::summarise(width = sum(width)) %>%
        as.data.frame()
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
}
















