#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create factRObject
#'
#' @description
#' factRObject stores all imported inputs, intermediate objects and results of
#' factR workflow. To construct a factRObject, you need a custom transcriptome
#' in GTF format, a reference annotation and a genome sequence. See examples
#' [NOT COMPLETE]
#'
#'
#' @slot ranges a list of GenomicRanges GTF objects
#' @slot genome genome sequence
#' @slot domains a dataFrame of predicted protein domains
#' @slot nmd a list of NMD-related objects
#' @slot misc  a list of miscellaneous information
#' @slot version version of factR this object was built under
#'
#' @return A factRObject
#' @name factRObject-class
#' @rdname factRObject-class
#' @exportClass factR
#'
setClass("factR",
         slots = c(
             ranges = "list",
             genome = "list",
             domains = "data.frame",
             nmd = "list",
             misc = "list",
             version = "character"
         )
)


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
#' gtf <- system.file("extdata", "sc_merged.gtf.gz", package = "factR")
#'
#' sample.factr <- CreatefactRObject(gtf, "vM25", Mmusculus)
CreatefactRObject <- function(gtf, annotation, genome){

    # catch missing args
    mandargs <- c("gtf", "annotation","genome")
    passed <- names(as.list(match.call())[-1])
    if (any(!mandargs %in% passed)) {
        rlang::abort(paste(
            "missing values for",
            paste(setdiff(mandargs, passed), collapse = ", ")
        ))
    }

    # check existence of GTF input and try to import
    if (file.exists(gtf) & tolower(tools::file_ext(gtf)) == "gtf"){
        rlang::abort(sprintf("File '%s' does not exist", gtf))
    } else if (!grepl("gtf", tolower(gtf))){
        rlang::abort(sprintf("File '%s' is not a GTF", gtf))
    } else {
        obj <- methods::new("factR")
        #rlang::inform("Importing custom GTF")
        obj@ranges$custom <- factR::importGTF(gtf)  # import GTF
    }

    # check for type of annotation and import
    data(genomes)
    if(is.character(annotation)){
        if (annotation %in% genomes$ID){
            #rlang::inform("Retrieving reference transcriptome")
            path_to_gtf_ftp <- genomes[genomes$ID == annotation,]$annotation
            temp.gtf <- tempfile()
            download.file(path_to_gtf_ftp, temp.gtf, quiet = TRUE)
            obj@ranges$ref <- factR::importGTF(temp.gtf)
        } else if (file.exists(annotation) & grepl(".gtf", tolower(annotation))){
            #rlang::inform("Importing annotation")
            obj@ranges$ref <- factR::importGTF(annotation)
        } else {
            rlang::abort(sprintf("Annotation '%s' does not exist",
                                 annotation))
        }
    } else if(is.object(annotation)){
        if (!is_gtf(annotation)){
            rlang::abort("Annotation object is not in GTF format")
        } else {
            obj@ranges$ref <- annotation
        }
    }

    if(is.character(genome)){
        if (genome %in% genomes$ID){
            #rlang::inform("Retrieving reference genome")
            options(timeout = 100000)
            path_to_gtf_ftp <- genomes[genomes$ID == genome,]$genome
            temp.fa <- tempfile()
            download.file(path_to_gtf_ftp, temp.fa, quiet = TRUE)
            obj@genome$seq <- factR::importFASTA(temp.fa)
        } else if (file.exists(genome) & grepl(".fa", tolower(genome))){
            #rlang::inform("Importing genome")
            obj@genome$seq <- factR::importFASTA(genome)
        } else {
            rlang::abort(sprintf("Genome '%s' does not exist",
                                 genome))
        }
    } else if(is.object(genome)){
        if (!class(genome) %in% c("BSgenome", "DNAStringSet")){
            rlang::abort("Genome object is not in the right format")
        } else {
            obj@genome$seq <- genome
        }
    }
    return(obj)
}




















