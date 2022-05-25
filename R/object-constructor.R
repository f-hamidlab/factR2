#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @param gtf path to custom GTF file.
#' @param annotation a reference transcriptome annotation. Input can be (1) a
#' GenomicRanges object containing transcriptome GTF information, (2) a path
#' to a reference GTF file on local computer or (3) a character value representing
#' one of the genome IDs that we support. Run \code{\link{ListGenomes}} for our
#' list of supported genomes
#' @param genome genome sequence. Input can be (1) a
#' BSgenome object containing genome sequence, (2) a path
#' to a reference GTF file on local computer or (3) a character value representing
#' one of the genome IDs that we support. Run \code{\link{ListGenomes}} for our
#' list of supported genomes
#'
#' @name factRObject-class
#' @rdname factRObject-class
#' @return factR object
#' @export
#'
#' @examples
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' gtf <- system.file("extdata", "sc_merged_sample.gtf.gz", package = "factR")
#'
#' sample.factr <- CreatefactRObject(gtf, "vM25", Mmusculus)
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
            rlang::inform("Retrieving custom reference")
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


















