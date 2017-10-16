\name{GeneSet-class}
\Rdversion{1.1}
\docType{class}
\alias{GeneSet-class}

\title{Class \code{"GeneSet"}}
\description{
Class for storing sets of genes and their corresponding metadata.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GeneSet")}.

These objects are used internally by the \code{chipenrich} package and users will not likely need to create these.
}

\section{Slots}{
  \describe{
    \item{\code{set.gene}:}{Object of class \code{"environment"}. Maps from geneset IDs to lists of Entrez gene IDs. }
    \item{\code{type}:}{Object of class \code{"character"}. The formal name for this collection of genesets. }
    \item{\code{set.name}:}{Object of class \code{"environment"}. Maps from geneset IDs to their descriptions/names. }
    \item{\code{all.genes}:}{Object of class \code{"character"}. A set of all genes present across every geneset. }
    \item{\code{organism}:}{Object of class \code{"character"}. Organism code for gene IDs. }
    \item{\code{dburl}:}{Object of class \code{"character"}. Web URL for this collection of genesets. }
  }
}

\author{
Ryan Welch \email{welchr@umich.edu}
}

\note{
Not typically accessed by the user - this is used internally by the \code{chipenrich} package.
}

\examples{
# Show information about class.
showClass("GeneSet")

# What is stored inside a geneset object?
data("geneset.metabolite.hsa");
str(geneset.metabolite.hsa);

# How are the mappings from geneset IDs to gene IDs stored?
ls.str(geneset.metabolite.hsa@set.gene);
}

\keyword{classes}
