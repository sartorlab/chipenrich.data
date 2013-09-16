\name{chipenrich.data-package}
\alias{chipenrich.data-package}
\alias{chipenrich.data}
\docType{package}

\title{
Datasets required for the chipenrich package
}

\description{
This package includes genesets (KEGG, GO, etc.), gene locus definitions, gene locus mappabilities, and other required data for the chipenrich package. It also includes example ChIP-seq peak datasets for use with chipenrich. 
}

\details{
Almost all of the datasets within this package, with the exception of the example ChIP-seq peak datasets, are not intended for direct use by the user. They are automatically used by the \code{chipenrich} package itself. This help file attempts to document each data source should a user wish to understand how each is computed, where it came from, etc. 

Genesets are exactly that - sets of genes - as defined in various databases such as KEGG, Gene Ontology, and others. 

Each geneset is stored per genome in the format "geneset.<name of geneset>.<genome build>", use \code{data(package="chipenrich.data")} for a full list. 

Each geneset and the description for how it was generated can be found in the "Geneset Database Information" vignette. Use \code{browseVignettes("chipenrich.data")} for a full list. 

Gene locus definitions define regions in which peaks are assigned to genes. If a peak falls into a particular gene locus, it is assigned to that gene. Within the package, locus definitions are stored per genome and per definition in the format "locusdef.<genome>.<definition code>". We support multiple definitions currently: 

\itemize{
	\item \code{nearest_tss}: assigns peaks to the gene with the closest transcription start site (TSS). Each gene locus stretches from the midpoint upstream between TSSs, and the midpoint downstream between TSSs. 
	\item \code{nearest_gene}: assigns peaks to the nearest gene. Each gene locus stretches from the upstream midpoint between genes and the downstream midpoint between genes. 
	\item \code{1kb}: each gene locus is defined as the region 1KB up- and down-stream of the TSS. 
	\item \code{5kb}: each gene locus is defined as the region 5KB up- and down-stream of the TSS. 
}

The mappability datasets stored within this package represent the average mappability per locus definition (above.) They are stored in the following format: "mappa.<genome build>.<locus definition code>.<read length>". Each dataset is simply an entrez ID, and the mappability score computed for the gene under the given locus definition and read length. For more information regarding how mappability is included, please see (Welch & Lee, 2012.) 

Example ChIP-seq peak datasets are also included, and can be used for running \code{chipenrich}. Each dataset is stored as "peaks_<name of protein/transcription factor>". Use \code{\link{example_datasets}} or \code{data(package="chipenrich.data")} for a list. 

}

\author{
Ryan Welch \email{welchr@umich.edu}
}

\references{
	Welch RP*, Lee C*, Smith RA, Imbriano P, Scott LJ, Sartor MA. 2012. ChIP-Enrich: Gene set enrichment testing for ChIP-Seq data (in preparation.)
}

\examples{
# Get a list of all ChIP-seq peak example datasets.
example_datasets();

# Inspect an example geneset.
data("geneset.kegg_pathway.hsa");
str(geneset.kegg_pathway.hsa);

# Get the mappings from geneset ID to genes for KEGG (loaded above)
kegg_to_genes = as.list(geneset.kegg_pathway.hsa@set.gene);

# Get the mappings from geneset ID to geneset descriptions for KEGG
kegg_to_names = as.list(geneset.kegg_pathway.hsa@set.name);

# List all genes within KEGG pathways (as Entrez IDs.)
kegg_genes = geneset.kegg_pathway.hsa@all.genes;

# Inspect an example locus definition.
data(locusdef.hg19.nearest_tss);
str(locusdef.hg19.nearest_tss);

# Get a data frame specifying each gene locus, its chromosome, and start/end positions.
ldef = locusdef.hg19.nearest_tss@dframe;

# Get a GRanges object containing each locus definition for fast lookups.
granges = locusdef.hg19.nearest_tss@granges;
}

