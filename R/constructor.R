#' Construct factR object class
#'
#' @description
#' Create a factR object from custom GTF transcriptomes
#'
#'
#'
#'
#' @param gtf Either path to custom transcriptome file (GTF) or a GenomicRanges
#' object containing GTF transcriptome data
#' @param reference Character value of the ID of genome used as reference. See
#' \link{listSupportedGenomes()} to print a list of supported genomes. Input can
#' also be name of species (Homo sapiens or Hsapiens or Human).
#' @param use_own_annotation Can be one of the following:
#' \itemize{
#'  \item{}{Path to local transcriptome file}
#'  \item{}{GenomicRanges object containing GTF transcriptome data}
#'  \item{}{AnnotationHub data ID [AHxxxxx]}
#'  \item{}{URL to GTF file}
#' }
#' @param use_own_genome Can be one of the following:
#' \itemize{
#'  \item{}{Path to local FASTA file}
#'  \item{}{Biostrings object containing full genome sequence}
#'  \item{}{AnnotationHub data ID [AHxxxxx]}
#'  \item{}{URL to FASTA file}
#' }
#' @param project_name Character value of the name of project
#' @param genome_build Character value of the genome build. Will be determined
#' automatically by default.
#' @param match_genes Boolean value as to whether genes in custom transcriptome
#' is to be matched to reference (Default: TRUE)
#' @param countData (Optional) Matrix of transcript-level counts data
#' @param sampleData (Optional) Dataframe containing sample metadata
#' @param verbose Boolean value as to whether messages should be printed (Default: TRUE)
#'
#' @return
#' factRObject class.
#'
#' @seealso \code{\link{factRObject-class}} \code{\link{factR-exp}}
#'
#' @export
#'
#' @examples
#' gtf <- system.file("extdata/sc_merged_sample.gtf.gz", package = "factR")
#' factR.object <- createfactRObject(gtf, "vM25")
#'
createfactRObject <- function(gtf, reference,
                              use_own_annotation = NULL,
                              use_own_genome = NULL,
                              project_name = "factRProject",
                              genome_build = "auto",
                              match_genes = TRUE,
                              countData = NULL,
                              sampleData = NULL,
                              psi = NULL,
                              verbose = TRUE){

    # check input arguments
    if(verbose){
        .msgheader("Checking inputs")
    }
    ## GTF input
    if(missing(gtf)){
        rlang::abort("Missing gtf input")
    }
    ## Annotation input
    reference.list <- .checkannotation(reference, use_own_annotation,
                                       use_own_genome, verbose)

    ## Count data input
    # if(!is.null(countData)){
    #     # check if sampleData and design are given
    #     if(is.null(sampleData)){
    #         rlang::abort("Counts data provided, but sample data missing")
    #     }
    # }


    # create new factRObject
    if(verbose){
        .msgheader("Checking factRObject")
    }
    obj <- methods::new("factR")
    obj@version <- factR2version
    obj@project <- project_name

    # import data
    if(verbose){
        .msgheader("Adding custom transcriptome")
    }
    obj@transcriptome <- .smartimport(gtf, ".gtf", "Custom transcriptome",verbose)
    GenomeInfoDb::seqlevels(obj@transcriptome) <- as.character(
        unique(GenomicRanges::seqnames(obj@transcriptome)))

    if(verbose){
        .msgheader("Adding annotation")
    }
    obj@reference$ranges <- .smartimport(reference.list[[1]], ".gtf", "Annotation",
                                         verbose)
    if(verbose){
        .msgheader("Adding genome sequence")
    }
    obj@reference$genome <- .smartimport(reference.list[[2]], ".fa", "Genome sequence",
                                         verbose)

    # set genome build
    if(genome_build == "auto"){
        obj@reference$build <- reference.list[[3]]
    } else {
        obj@reference$build <- genome_build
    }


    # run object prepration function
    obj <- .prepfactR(obj, match_genes, verbose)

    # add and prep counts data if given
    if(!is.null(countData)){
        obj <- addTxCounts(obj, countData, sampleData)
    }
    .msginfo("factRobject created!")
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
            build <- "custom"
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
            suppressPackageStartupMessages(require(selected_genome$genome.pri, character.only = T,
                    quietly = !verbose))
            genome <- get(selected_genome$genome.pri)
        } else {
            # prompt to download BSgenome
            .msgsubinfo(stringr::str_glue(
                "The BSgenome object {selected_genome$genome.pri} is available
                for this annotation. Do you wish to download? [Y/N]"))
            resp <- stringr::str_to_upper(readline(prompt=""))
            if(resp=="Y"){
                suppressMessages(BiocManager::install(selected_genome$genome.pri,
                                     quiet=TRUE,
                                     ask=FALSE))
                genome <- get(selected_genome$genome.pri)

            } else{
                genome <- selected_genome$genome.sec
            }


        }
        build <- selected_genome$build
    }
    return(list(annotation, genome, build))

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
            .msgsubinfo(sprintf("Using %s object", datatype))
        }
        return(input)
    }



    ## try local file
    else if (file.exists(input) & grepl(format, tolower(input))){
        if(verbose){
            .msgsubinfo("Importing from local directory")
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
            .msgsubinfo("Importing from AnnotationHub")

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
            .msgsubinfo("Importing from URL")
        }
        if(format == ".gtf"){
            tmpfile <- tempfile()
            download.file(input, tmpfile, quiet = TRUE)
            return(factR::importGTF(tmpfile))
        } else if(format == ".fa"){
            tmpfile <- tempfile()
            download.file(input, tmpfile, quiet = TRUE)
            return(factR::importFASTA(tmpfile))
        }

    } else {
        rlang::abort(sprintf("Input not recognized or file does not exists"))
    }
}



.prepfactR <-  function(object, matchgenes, verbose = FALSE) {

    # check chromosome overlap
    `your custom GTF` <- object@transcriptome
    `reference annotation` <- object@reference$ranges
    result <- factR::has_consistentSeqlevels(`your custom GTF` ,
                                             `reference annotation`)


    #  match chromosomes
    if(verbose){
        .msgheader("Matching chromosome names")
    }
    object@transcriptome <- suppressWarnings(
        factR::matchChromosomes(object@transcriptome,
                                                    object@reference$genome))
    object@reference$ranges <- suppressWarnings(
        factR::matchChromosomes(object@reference$ranges,
                                object@reference$ranges))

    # match gene ID if requested
    if(matchgenes){
        ## try find variables that contain gene ids
        if(verbose){
            .msgheader("Matching gene names")
        }
        if("ref_gene_id" %in% colnames(as.data.frame(object@transcriptome))){
            if(verbose){
                object@transcriptome <- suppressWarnings(
                    factR::matchGeneInfo(object@transcriptome,
                                         object@reference$ranges,
                                         primary_gene_id = "gene_id",
                                         secondary_gene_id = "ref_gene_id"))
} else {
                object@transcriptome <- suppressMessages(suppressWarnings(
                    factR::matchGeneInfo(object@transcriptome,
                                         object@reference$ranges,
                                         primary_gene_id = "gene_id",
                                         secondary_gene_id = "ref_gene_id")))
            }
        } else {
            if(verbose){
                object@transcriptome <- suppressWarnings(
                    factR::matchGeneInfo(object@transcriptome,
                                             object@reference$ranges,
                                             primary_gene_id = "gene_id"))
            } else {
                object@transcriptome <- suppressWarnings(suppressMessages(
                    factR::matchGeneInfo(object@transcriptome,
                                         object@reference$ranges,
                                         primary_gene_id = "gene_id")))
            }
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
        .msgheader("Creating factRset objects")
    }
    if(verbose){
        .msgsubinfo("Adding gene information")
    }
    object@active.set <- "AS"
    object@sets$gene <- methods::new("factRset")
    object@sets$gene@rowData <- as.data.frame(object@transcriptome) %>%
        dplyr::filter(type %in% "gene") %>%
        dplyr::select(gene_id, gene_name, strand, width, match_level) %>%
        dplyr::distinct()
    rownames(object@sets$gene@rowData) <- object@sets$gene@rowData$gene_id
    object@sets$gene@counts <- as.matrix(data.frame(row.names =  rownames(object[["gene"]])))
    object@sets$gene@data <- as.matrix(data.frame(row.names =  rownames(object[["gene"]])))

    if(verbose){
        .msgsubinfo("Adding transcript information")
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
        .msgsubinfo("Adding alternative splicing information")
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
        .msgheader("Annotating novel transcripts")
    }
    gtf <- granges(object, set = "transcript")
    newtxs <- suppressMessages(suppressWarnings(
        factR::subsetNewTranscripts(gtf,
                                    object@reference$ranges,
                                    refine.by = "intron")))
    object <- addMeta(object,
                      meta="transcript",
                      novel = ifelse(transcript_id %in% newtxs$transcript_id,
                                     "yes","no"),
                      cds = "no",
                      nmd = "no")

    # annotate new events
    if(verbose){
        .msgheader("Annotating novel AS events")
    }
    ref.gtf <- object@reference$ranges
    ref.AS <- .runAS(ref.gtf[ref.gtf$type == "exon"])

    AS.id <- ase(object, show_more = TRUE) %>%
        dplyr::mutate(id = paste0(coord,gene_id,strand,AStype)) %>%
        dplyr::mutate(id = ifelse(AStype == "AF" & strand == "+",
                                  paste0(stringr::str_replace(
                                      coord, ":[0-9]+-",":"),
                                      gene_id,strand,AStype),
                                  id)) %>%
        dplyr::mutate(id = ifelse(AStype == "AF" & strand == "-",
                                  paste0(stringr::str_replace(
                                      coord, "-[0-9]+$",""),
                                      gene_id,strand,AStype),
                                  id)) %>%
        dplyr::mutate(id = ifelse(AStype == "AL" & strand == "+",
                                  paste0(stringr::str_replace(
                                      coord, "-[0-9]+$",""),
                                      gene_id,strand,AStype),
                                  id))  %>%
        dplyr::mutate(id = ifelse(AStype == "AL" & strand == "-",
                                  paste0(stringr::str_replace(
                                      coord, ":[0-9]+-",":"),
                                      gene_id,strand,AStype),
                                  id)) %>%
        dplyr::pull(id)
    ref.AS.id <- ref.AS %>%
        as.data.frame() %>%
        dplyr::mutate(coord = paste0(seqnames, ":", start, "-", end)) %>%
        dplyr::mutate(id = paste0(coord,gene_id,strand,AStype)) %>%
        dplyr::mutate(id = ifelse(AStype == "AF" & strand == "+",
                                  paste0(stringr::str_replace(
                                      coord, ":[0-9]+-",":"),
                                      gene_id,strand,AStype),
                                  id)) %>%
        dplyr::mutate(id = ifelse(AStype == "AF" & strand == "-",
                                  paste0(stringr::str_replace(
                                      coord, "-[0-9]+$",""),
                                      gene_id,strand,AStype),
                                  id)) %>%
        dplyr::mutate(id = ifelse(AStype == "AL" & strand == "+",
                                  paste0(stringr::str_replace(
                                      coord, "-[0-9]+$",""),
                                      gene_id,strand,AStype),
                                  id))  %>%
        dplyr::mutate(id = ifelse(AStype == "AL" & strand == "-",
                                  paste0(stringr::str_replace(
                                      coord, ":[0-9]+-",":"),
                                      gene_id,strand,AStype),
                                  id)) %>%
        dplyr::pull(id)

    object@sets$AS@rowData$novel <- ifelse(!AS.id %in% ref.AS.id, "yes", "no")
    object
}
















