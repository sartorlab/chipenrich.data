\name{mappa.mm9.1kb.100mer}

\alias{mappa.mm9.1kb.100mer}

\docType{data}

\title{
mappa.mm9.1kb.100mer
}

\description{
Gene locus mappability data, calculated for: 

\itemize{
  \item Build: mm9
  \item Locus definition: 1kb
  \item K-Mer Reads: 100mer   
}

These objects are used internally by the \code{chipenrich} package. 
}

\usage{data(mappa.mm9.1kb.100mer)}

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
