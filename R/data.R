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

