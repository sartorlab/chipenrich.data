\name{mappa.hg19.1kb.24mer}

\alias{mappa.hg19.1kb.24mer}

\docType{data}

\title{
mappa.hg19.1kb.24mer
}

\description{
Gene locus mappability data, calculated for: 

\itemize{
  \item Build: hg19
  \item Locus definition: 1kb
  \item K-Mer Reads: 24mer   
}

These objects are used internally by the \code{chipenrich} package. 
}

\usage{data(mappa.hg19.1kb.24mer)}

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