% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructor.R
\name{createfactRObject}
\alias{createfactRObject}
\title{Construct factR object class}
\usage{
createfactRObject(
  gtf,
  reference,
  use_own_annotation = NULL,
  use_own_genome = NULL,
  project_name = "factRProject",
  genome_build = "auto",
  match_genes = TRUE,
  countData = NULL,
  sampleData = NULL,
  psi = NULL,
  verbose = TRUE
)
}
\arguments{
\item{gtf}{Either path to custom transcriptome file (GTF) or a GenomicRanges
object containing GTF transcriptome data}

\item{reference}{Character value of the ID of genome used as reference. See
\link{listSupportedGenomes()} to print a list of supported genomes. Input can
also be name of species (Homo sapiens or Hsapiens or Human).}

\item{use_own_annotation}{Can be one of the following:
\itemize{
 \item{}{Path to local transcriptome file}
 \item{}{GenomicRanges object containing GTF transcriptome data}
 \item{}{AnnotationHub data ID [AHxxxxx]}
 \item{}{URL to GTF file}
}}

\item{use_own_genome}{Can be one of the following:
\itemize{
 \item{}{Path to local FASTA file}
 \item{}{Biostrings object containing full genome sequence}
 \item{}{AnnotationHub data ID [AHxxxxx]}
 \item{}{URL to FASTA file}
}}

\item{project_name}{Character value of the name of project}

\item{genome_build}{Character value of the genome build. Will be determined
automatically by default.}

\item{match_genes}{Boolean value as to whether genes in custom transcriptome
is to be matched to reference (Default: TRUE)}

\item{countData}{(Optional) Matrix of counts data}

\item{sampleData}{Dataframe containing sample metadata}

\item{psi}{(Optional) Matrix of splicing PSI data}

\item{verbose}{Boolean value as to whether messages should be printed (Default: TRUE)}
}
\value{
factR object class
}
\description{
Create a factR object from custom GTF transcriptome
}
\examples{
gtf <- system.file("extdata/sc_merged_sample.gtf.gz", package = "factR")
factR.object <- createfactRObject(gtf, "mm10")

}