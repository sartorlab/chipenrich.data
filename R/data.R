#' tss.hg19_refseq TSS locations
#'
#' A \code{GRanges} with all the TSSs for hg19_refseq. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: UCSC_refGene_table_02072018 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"tss.hg19_refseq"

