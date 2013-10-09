\name{mappa.mm9.5kb.40mer}

\alias{mappa.mm9.5kb.40mer}

\docType{data}

\title{
mappa.mm9.5kb.40mer
}

\description{
Gene locus mappability data, calculated for: 

\itemize{
  \item Build: mm9
  \item Locus definition: 5kb
  \item K-Mer Reads: 40mer   
}

These objects are used internally by the \code{chipenrich} package. 
}

\usage{data(mappa.mm9.5kb.40mer)}

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