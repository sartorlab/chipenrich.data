\name{LocusDefinition-class}
\Rdversion{1.1}
\docType{class}
\alias{LocusDefinition-class}

\title{Class \code{"LocusDefinition"}}
\description{
A storage class representing gene locus definitions and their corresponding metadata.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("LocusDefinition")}.
These objects are used internally by the \code{chipenrich} package and users will not likely need to create these.
}

\section{Slots}{
  \describe{
    \item{\code{dframe}:}{Object of class \code{"data.frame"}. Each row represents a locus for a particular geneid. }
    \item{\code{granges}:}{Object of class \code{"GenomicRanges"}. Locus definitions stored as a GenomicRanges object. }
    \item{\code{genome.build}:}{Object of class \code{"character"}. Genome build these definitions were generated from. }
    \item{\code{organism}:}{Object of class \code{"character"}. Organism code.  }
  }
}

\author{
Ryan Welch \email{welchr@umich.edu}
}

\note{
Not typically accessed by the user - this is used internally by the \code{chipenrich} package.
}

\examples{
# Show info about the class.
showClass("LocusDefinition");

# Example of what a locus definition object looks like.
data(locusdef.mm9.nearest_tss)
lapply(attributes(locusdef.mm9.nearest_tss), head)
}

\keyword{classes}
