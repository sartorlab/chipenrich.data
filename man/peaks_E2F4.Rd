\name{peaks_E2F4}
\docType{data}
\alias{peaks_E2F4}

\title{
  ChIP-seq Peaks for the E2F4 Transcription Factor
}

\description{
	A dataset containing the binding locations (peaks) of the transcription factor E2F4 called from a ChIP-seq experiment
}

\usage{
peaks_E2F4
}

\format{
  A data frame containing 16,245 peak binding locations. Each row is a peak location, and the 3 variables/columns are chromosome, start peak coordinate, and end peak coordinate.
	Peak coordinates are in hg19 (UCSC) coordinates.
}

\source{
  The data and information regarding the experiment can be found in the following publication:

	Lee, B. K., A. A. Bhinge, et al. (2011). "Wide-ranging functions of E2F4 in transcriptional activation and repression revealed by genome-wide analysis." Nucleic Acids Res 39(9): 3558-3573.
}

\examples{
# Load E2F4 peak data. 
data(peaks_E2F4)

# Print the first 10 peaks in the dataset.
print(head(peaks_E2F4))
}
\keyword{datasets}
