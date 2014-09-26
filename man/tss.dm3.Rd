\name{tss.dm3}
\alias{tss.dm3}
\docType{data}

\title{
tss.dm3
}

\description{
Transcription start site information for dm3. 

This is used internally by the \code{chipenrich} package. 
}

\usage{data(tss.dm3)}

\format{
  A list containing: 
  \describe{
    \item{\code{data}}{A data frame containing TSSs for each gene.}
    \item{\code{ranges}}{A list mapping chromosomes to IRanges objects.}
    \item{\code{granges}}{TSSs stored as a GenomicRanges object, including strand info.}
  }
}

\author{
Chee Lee \email{cheelee@umich.edu}
}

\seealso{
\code{\link{chipenrich.data}}
}

\keyword{datasets}
