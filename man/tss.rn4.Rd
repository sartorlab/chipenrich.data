\name{tss.rn4}
\alias{tss.rn4}
\docType{data}

\title{
tss.rn4
}

\description{
Transcription start site information for rn4. 

This is used internally by the \code{chipenrich} package. 
}

\usage{data(tss.rn4)}

\format{
  A list containing: 
  \describe{
    \item{\code{data}}{A data frame containing TSSs for each gene.}
    \item{\code{ranges}}{A list mapping chromosomes to IRanges objects.}
    \item{\code{granges}}{TSSs stored as a GenomicRanges object, including strand info.}
  }
}

\author{
Ryan Welch \email{welchr@umich.edu}
}

\seealso{
\code{\link{chipenrich.data}}
}

\keyword{datasets}
