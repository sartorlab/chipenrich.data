% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genesets_GO_doc.R
\docType{data}
\name{geneset.GOCC.dre}
\alias{geneset.GOCC.dre}
\title{geneset.GOCC.dre genesets for Danio rerio}
\format{A \code{GeneSet} object with the following slots:
\describe{
    \item{type}{A \code{character} indicating the type of genesets, e.g. GOBP.}
    \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
    \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
    \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. GO IDs for GO terms), and whose elements are \code{character} vectors of Entrez Gene IDs.}
    \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
    \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
}}
\source{
org.Dr.eg.db_3.4.2 and GO.db_3.4.2
}
\usage{
geneset.GOCC.dre
}
\description{
Gene Ontology Cellular Component (GOCC) genesets for Danio rerio. All genesets are required to have >= 10 Entrez IDs.
Built on Mon Oct 16 16:25:08 2017.
}
\keyword{datasets}
