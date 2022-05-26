#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
#' @name factRObject-class
#' @rdname factRObject-class
#' @return
#' @export
#'
#' @examples
#' # Create factRObject using sample GTF and supported reference
#' gtf <- system.file("extdata", "sc_merged_sample.gtf.gz", package = "factR")
#' obj <- CreatefactRObject(gtf, "vM25")
#'
#' # Create factRObject using custom reference
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' obj <- CreatefactRObject(gtf, use_own_annotation = "AH49547", use_own_genome = Mmusculus)
CreatefactRObject <- function(gtf, reference,
                              use_own_annotation = NULL,
                              use_own_genome = NULL,
                              verbose = FALSE){

    # catch missing args
    mandargs <- c("gtf")
    passed <- names(as.list(match.call())[-1])
    if (any(!mandargs %in% passed)) {
        rlang::abort(paste(
            "missing values for",
            paste(setdiff(mandargs, passed), collapse = ", ")
        ))
    }
    if(missing(reference)){
        empty_custom_ref <- c(is.null(use_own_annotation), is.null(use_own_genome))
        if(all(empty_custom_ref)){
            rlang::abort("missing value for reference")
        } else if(any(empty_custom_ref)){
            rlang::abort(sprintf("Using own reference, but value for '%s' not provided",
                c("use_own_annotation", "use_own_genome")[empty_custom_ref]))
        }
    }

    # check existence of GTF input and try to import
    if (!file.exists(gtf)){
        rlang::abort(sprintf("File '%s' does not exist", gtf))
    } else if (!grepl("gtf", tolower(gtf))){
        rlang::abort(sprintf("File '%s' is not a GTF", gtf))
    } else {
        obj <- methods::new("factR")
        obj@version <- factR2version
        if(verbose){
            rlang::inform("Importing custom GTF")
        }
        obj@custom$ranges <- factR::importGTF(gtf)  # import GTF
    }


    # check for type of annotation and import
    data(genomes)
    # import supported reference
    if(!missing(reference)){

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

        # get annotation
        if(verbose){
            rlang::inform("Retrieving supported reference")
        }
        custom.gtf <- .smartimport(selected_genome$annotation, ".gtf", verbose)
        if(is_gtf(custom.gtf)){
            obj@reference$ranges <- custom.gtf
        } else {
            rlang::abort("Annotation object is not a GTF")
        }

        # get genome
        ## test if BSgenome object is installed
        if(selected_genome$genome.pri %in% BSgenome::installed.genomes()){
            require(selected_genome$genome.pri, character.only = T,
                    quietly = !verbose)
            obj@reference$genome <- get(selected_genome$genome.pri)
        } else {
            obj@reference$genome <- .smartimport(selected_genome$genome.sec,
                                                 ".fa",
                                                 verbose)
        }

        obj@reference$id <- selected_genome$ID
        obj@reference$species <- selected_genome$species
        obj@reference$db <- selected_genome$database
    # import custom reference
    } else {
        # get custom annotation
        if(verbose){
            rlang::inform("Retrieving user-specified reference")
        }

        custom.gtf <- .smartimport(use_own_annotation, ".gtf", verbose)
        if(is_gtf(custom.gtf)){
            obj@reference$ranges <- custom.gtf
        } else {
            rlang::abort("Annotation object is not a GTF")
        }

        # get custom genome
        obj@reference$genome <- .smartimport(use_own_genome, ".fa",
                                             verbose)
        obj@reference$id <- "custom"
        obj@reference$species <- "custom"
        obj@reference$db <- "custom"
    }

    # match chromosomes and gene ID
    if(verbose){
        rlang::inform("Matching chromosome names")
    }
    obj@custom$ranges <- factR::matchChromosomes(obj@custom$ranges,
                                                 obj@reference$genome)
    obj@reference$ranges <- suppressWarnings(factR::matchChromosomes(obj@reference$ranges,
                                                 obj@reference$ranges))
    ## try find variables that contain gene ids
    if(verbose){
        rlang::inform("Matching gene names")
    }
    custom.df <- as.data.frame(obj@custom$ranges)
    potential_id_vars <- apply(custom.df, 2, function(x) any(grepl("ENS",x)))
    potential_id <- names(potential_id_vars)[potential_id_vars]
    potential_id <- potential_id[-which("transcript_id" %in% potential_id)]
    if(length(potential_id > 1)){
        if(verbose){
            obj@custom$ranges <- factR::matchGeneInfo(obj@custom$ranges,
                                                       obj@reference$ranges,
                                                       primary_gene_id = "gene_id",
                                                       secondary_gene_id = potential_id)
        } else {
            obj@custom$ranges <- suppressMessages(
                factR::matchGeneInfo(obj@custom$ranges,
                                     obj@reference$ranges,
                                     primary_gene_id = "gene_id",
                                     secondary_gene_id = potential_id))
        }

    } else {
        if(verbose){
            obj@custom$ranges <- factR::matchGeneInfo(obj@custom$ranges,
                                                      obj@reference$ranges,
                                                      primary_gene_id = "gene_id")
        } else {
            obj@custom$ranges <- suppressMessages(
                factR::matchGeneInfo(obj@custom$ranges,
                                     obj@reference$ranges,
                                     primary_gene_id = "gene_id"))
        }
    }

    # create genetxs dataframe
    if(verbose){
        rlang::inform("Creating a list of identified transcripts")
    }
    obj@custom$genetxs <- as.data.frame(obj@custom$ranges) %>%
        dplyr::select(gene_id, gene_name, transcript_id, match_level) %>%
        dplyr::distinct()

    # annotate new transcripts
    newtxs <- suppressMessages(factR::subsetNewTranscripts(obj@custom$ranges,
                                                           obj@reference$ranges,
                                                           refine.by = "intron"))
    obj@custom$genetxs$novel <- ifelse(obj@custom$genetxs$transcript_id %in% newtxs$transcript_id,
                                       "yes", "no")
    obj@custom$genetxs$cds <- "no"
    obj@custom$genetxs$nmd <- "no"

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


















