\name{mappa.hg19.nearest_tss.75mer}

\alias{mappa.hg19.nearest_tss.75mer}

\docType{data}

\title{
mappa.hg19.nearest_tss.75mer
}

\description{
Gene locus mappability data, calculated for: 

\itemize{
  \item Build: hg19
  \item Locus definition: nearest_tss
  \item K-Mer Reads: 75mer   
}

These objects are used internally by the \code{chipenrich} package. 
}

\usage{data(mappa.hg19.nearest_tss.75mer)}

\format{
  A data frame containing: 
  \describe{
    \item{\code{geneid}}{Entrez Gene IDs}
    \item{\code{mappa}}{Gene locus mappability}
  }
}

\author{
Ryan Welch \email{welchr@umich.edu}
}

\seealso{
For more information about gene locus definitions: \code{\link{chipenrich.data}}
}

\keyword{datasets}
