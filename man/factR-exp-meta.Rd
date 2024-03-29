% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/factRObject-class.R
\name{factR-exp-meta}
\alias{factR-exp-meta}
\alias{samples,factR-method}
\alias{colNames,factR-method}
\alias{colNames<-,factR-method}
\alias{ident,factR-method}
\alias{counts,factR-method}
\title{Display factR2 samples and expression data}
\usage{
\S4method{rowNames}{factR}(object)

\S4method{samples}{factR}(object)

\S4method{colNames}{factR}(object)

\S4method{colNames}{factR}(object) <- value

\S4method{ident}{factR}(object)

\S4method{counts}{factR}(object, ..., set = NULL, slot = "data")
}
\arguments{
\item{object}{factRObject}

\item{value}{New sample names. Length should match number of columns}

\item{...}{One or more features to display. Can be the following:
\itemize{
 \item{gene_id: }{ID of gene to plot}
 \item{gene_name: }{Name of gene to plot}
 \item{transcript_id: }{ID of transcript to plot}
}}

\item{set}{Set metadata to display. Can be "gene", "transcript" or "AS".}

\item{slot}{Data slot to display. Can be "counts" or "data" (default)}
}
\value{
\itemize{
 \item{`counts`: }{Count matrix}
 \item{All other functions: }{Character vector of values}
}
}
\description{
If expression data from multiple samples is provided,
factRObject-class will store the expression values at the gene, transcript and
alternative exon levels as well as samples metadata.
The functions below display contents pertaining to expression and samples data:
}
\details{
See \code{\link{factR-exp}} for example codes
}
\seealso{
\code{\link{factR-exp}}
}
