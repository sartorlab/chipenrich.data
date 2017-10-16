\name{peaks_H3K4me3_GM12878}
\docType{data}
\alias{peaks_H3K4me3_GM12878}

\title{
  ChIP-seq Peaks for the Histone Modification H3K4me3 in GM12878
}

\description{
	A dataset containing the binding locations (peaks) of the histone modification H3K4me3 called from a ChIP-seq experiment in the GM12878 cell line
}

\usage{
peaks_H3K4me3_GM12878
}

\format{
  A data frame containing 57,476 peak binding locations. Each row is a peak location, and the 3 variables/columns are chromosome, start peak coordinate, and end peak coordinate.
	Peak coordinates are in hg19 (UCSC) coordinates.
}

\source{
  The data and information regarding the experiment can be found at the following ENCODE URL:

	http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneGm12878H3k04me3StdPkV2.broadPeak.gz
}

\examples{
# Load H3K4me3 in GM12878 peak data. 
data(peaks_H3K4me3_GM12878)

# Print the first 10 peaks in the dataset.
print(head(peaks_H3K4me3_GM12878))
}
\keyword{datasets}
