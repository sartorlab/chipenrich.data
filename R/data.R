#' tss.hg19 TSS locations
#'
#' A \code{GRanges} with all the TSSs for hg19. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"tss.hg19"

#' tss.hg38 TSS locations
#'
#' A \code{GRanges} with all the TSSs for hg38. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"tss.hg38"

#' tss.mm9 TSS locations
#'
#' A \code{GRanges} with all the TSSs for mm9. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"tss.mm9"

#' tss.mm10 TSS locations
#'
#' A \code{GRanges} with all the TSSs for mm10. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"tss.mm10"

#' tss.rn4 TSS locations
#'
#' A \code{GRanges} with all the TSSs for rn4. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.5.0.
"tss.rn4"

#' tss.rn5 TSS locations
#'
#' A \code{GRanges} with all the TSSs for rn5. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.2 and org.Rn.eg.db_3.5.0.
"tss.rn5"

#' tss.rn6 TSS locations
#'
#' A \code{GRanges} with all the TSSs for rn6. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.1 and org.Rn.eg.db_3.5.0.
"tss.rn6"

#' tss.dm3 TSS locations
#'
#' A \code{GRanges} with all the TSSs for dm3. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.5.0.
"tss.dm3"

#' tss.dm6 TSS locations
#'
#' A \code{GRanges} with all the TSSs for dm6. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.4.1 and org.Dm.eg.db_3.5.0.
"tss.dm6"

#' tss.danRer10 TSS locations
#'
#' A \code{GRanges} with all the TSSs for danRer10. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Drerio.UCSC.danRer10.refGene_3.4.2 and org.Dr.eg.db_3.5.0.
"tss.danRer10"

#' Enhancer locations
#'
#' A \code{GRanges} with all the enhancer locations for hg19. The locations were found using a combination of DNAse data and from Thurman et al (PMID: 22955617)
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
"enhancer.dnase_thurman.0"

#' Gene-Enhancer descriptives
#'
#' A data frame with gene-level descriptions of enhancer properties using enhancers.dnase_thurman.0. Used in the adjustment of proximity test to enhancers.
#' \describe{
#'     \item{gene_id}{The Entrez ID for the a gene}
#'     \item{avg_denh_emp}{The empirical average distance to an enhancer from 90 ENCODE ChIP-seq datasets. This is used as the adjustment.}
#'     \item{num_enh}{The number of enhancers assigned to the gene, defined by closest gene TSS}
#'     \item{avgdenh}{The theoretical average distance to an enhancer assuming every base pair on the genome is equally likely to have a peak binding.}
#' }
#'
"gene.enh.desc"

#' DTSS Spline adjustment
#'
#' A mgcv::gam object on a combined data of 90 ENCODE ChIP-seq datasets that modeled the relationship between a gene's locus length the distance from a peak to the gene's transcription start site, using a cubic spline. This is used to adjust for the proximity to TSSes test.
#'
"spline.log_dtss.90ENCODE"
