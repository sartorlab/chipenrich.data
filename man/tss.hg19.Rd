\name{tss.hg19}
\alias{tss.hg19}
\docType{data}

\title{
tss.hg19
}

\description{
Transcription start site information for hg19. 

This is used internally by the \code{chipenrich} package. 
}

\usage{data(tss.hg19)}

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
