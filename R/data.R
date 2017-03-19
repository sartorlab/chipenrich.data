#' tss.hg19 TSS locations
#'
#' A \code{GRanges} with all the TSSs for hg19. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"tss.hg19"

#' locusdef.hg19.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' Built on Sat Mar 18 12:52:09 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.nearest_tss"

#' locusdef.hg19.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' Built on Sat Mar 18 12:52:10 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.nearest_gene"

#' locusdef.hg19.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' Built on Sat Mar 18 12:52:12 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.exon"

#' locusdef.hg19.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' Built on Sat Mar 18 12:52:14 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.intron"

#' locusdef.hg19.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 12:52:15 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.1kb"

#' locusdef.hg19.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 12:52:15 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.1kb_outside_upstream"

#' locusdef.hg19.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 12:52:16 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.1kb_outside"

#' locusdef.hg19.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 12:52:17 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.5kb"

#' locusdef.hg19.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 12:52:17 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.5kb_outside_upstream"

#' locusdef.hg19.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 12:52:18 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.5kb_outside"

#' locusdef.hg19.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 12:52:19 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.10kb"

#' locusdef.hg19.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 12:52:19 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.10kb_outside_upstream"

#' locusdef.hg19.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 12:52:20 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.10kb_outside"

#' tss.hg38 TSS locations
#'
#' A \code{GRanges} with all the TSSs for hg38. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"tss.hg38"

#' locusdef.hg38.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' Built on Sat Mar 18 12:56:26 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.nearest_tss"

#' locusdef.hg38.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' Built on Sat Mar 18 12:56:27 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.nearest_gene"

#' locusdef.hg38.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' Built on Sat Mar 18 12:56:29 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.exon"

#' locusdef.hg38.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' Built on Sat Mar 18 12:56:31 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.intron"

#' locusdef.hg38.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 12:56:32 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.1kb"

#' locusdef.hg38.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 12:56:33 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.1kb_outside_upstream"

#' locusdef.hg38.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 12:56:34 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.1kb_outside"

#' locusdef.hg38.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 12:56:34 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.5kb"

#' locusdef.hg38.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 12:56:35 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.5kb_outside_upstream"

#' locusdef.hg38.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 12:56:36 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.5kb_outside"

#' locusdef.hg38.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 12:56:36 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.10kb"

#' locusdef.hg38.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 12:56:36 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.10kb_outside_upstream"

#' locusdef.hg38.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 12:56:37 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.10kb_outside"

#' tss.mm9 TSS locations
#'
#' A \code{GRanges} with all the TSSs for mm9. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"tss.mm9"

#' locusdef.mm9.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' Built on Sat Mar 18 12:57:30 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.nearest_tss"

#' locusdef.mm9.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' Built on Sat Mar 18 12:57:30 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.nearest_gene"

#' locusdef.mm9.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' Built on Sat Mar 18 12:57:32 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.exon"

#' locusdef.mm9.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' Built on Sat Mar 18 12:57:34 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.intron"

#' locusdef.mm9.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 12:57:35 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.1kb"

#' locusdef.mm9.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 12:57:35 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.1kb_outside_upstream"

#' locusdef.mm9.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 12:57:36 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.1kb_outside"

#' locusdef.mm9.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 12:57:36 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.5kb"

#' locusdef.mm9.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 12:57:37 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.5kb_outside_upstream"

#' locusdef.mm9.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 12:57:37 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.5kb_outside"

#' locusdef.mm9.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 12:57:38 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.10kb"

#' locusdef.mm9.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 12:57:38 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.10kb_outside_upstream"

#' locusdef.mm9.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 12:57:38 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.10kb_outside"

#' tss.mm10 TSS locations
#'
#' A \code{GRanges} with all the TSSs for mm10. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"tss.mm10"

#' locusdef.mm10.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' Built on Sat Mar 18 12:58:37 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.nearest_tss"

#' locusdef.mm10.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' Built on Sat Mar 18 12:58:37 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.nearest_gene"

#' locusdef.mm10.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' Built on Sat Mar 18 12:58:40 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.exon"

#' locusdef.mm10.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' Built on Sat Mar 18 12:58:41 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.intron"

#' locusdef.mm10.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 12:58:42 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.1kb"

#' locusdef.mm10.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 12:58:42 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.1kb_outside_upstream"

#' locusdef.mm10.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 13:06:19 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.1kb_outside"

#' locusdef.mm10.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 13:06:20 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.5kb"

#' locusdef.mm10.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 13:06:20 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.5kb_outside_upstream"

#' locusdef.mm10.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 13:06:21 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.5kb_outside"

#' locusdef.mm10.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 13:06:22 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.10kb"

#' locusdef.mm10.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 13:06:22 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.10kb_outside_upstream"

#' locusdef.mm10.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 13:06:23 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.4.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.10kb_outside"

#' tss.rn4 TSS locations
#'
#' A \code{GRanges} with all the TSSs for rn4. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.4.0.
"tss.rn4"

#' locusdef.rn4.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:06:44 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.4.0.
"locusdef.rn4.nearest_tss"

#' locusdef.rn4.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:06:44 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.4.0.
"locusdef.rn4.nearest_gene"

#' locusdef.rn4.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:06:46 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.4.0.
"locusdef.rn4.exon"

#' locusdef.rn4.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:06:47 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.4.0.
"locusdef.rn4.intron"

#' locusdef.rn4.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:06:47 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.4.0.
"locusdef.rn4.1kb"

#' locusdef.rn4.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:06:48 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.4.0.
"locusdef.rn4.1kb_outside_upstream"

#' locusdef.rn4.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:06:48 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.4.0.
"locusdef.rn4.1kb_outside"

#' locusdef.rn4.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:06:48 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.4.0.
"locusdef.rn4.5kb"

#' locusdef.rn4.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:06:49 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.4.0.
"locusdef.rn4.5kb_outside_upstream"

#' locusdef.rn4.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:06:49 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.4.0.
"locusdef.rn4.5kb_outside"

#' locusdef.rn4.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:06:49 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.4.0.
"locusdef.rn4.10kb"

#' locusdef.rn4.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:06:50 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.4.0.
"locusdef.rn4.10kb_outside_upstream"

#' locusdef.rn4.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:06:50 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.4.0.
"locusdef.rn4.10kb_outside"

#' tss.rn5 TSS locations
#'
#' A \code{GRanges} with all the TSSs for rn5. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"tss.rn5"

#' locusdef.rn5.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' Built on Sat Mar 18 13:07:06 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn5.nearest_tss"

#' locusdef.rn5.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' Built on Sat Mar 18 13:07:06 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn5.nearest_gene"

#' locusdef.rn5.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' Built on Sat Mar 18 13:07:08 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn5.exon"

#' locusdef.rn5.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' Built on Sat Mar 18 13:07:09 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn5.intron"

#' locusdef.rn5.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 13:07:10 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn5.1kb"

#' locusdef.rn5.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 13:07:10 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn5.1kb_outside_upstream"

#' locusdef.rn5.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 13:07:11 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn5.1kb_outside"

#' locusdef.rn5.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 13:07:11 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn5.5kb"

#' locusdef.rn5.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 13:07:11 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn5.5kb_outside_upstream"

#' locusdef.rn5.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 13:07:12 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn5.5kb_outside"

#' locusdef.rn5.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 13:07:12 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn5.10kb"

#' locusdef.rn5.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 13:07:12 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn5.10kb_outside_upstream"

#' locusdef.rn5.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 13:07:13 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn5.10kb_outside"

#' tss.rn6 TSS locations
#'
#' A \code{GRanges} with all the TSSs for rn6. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"tss.rn6"

#' locusdef.rn6.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' Built on Sat Mar 18 13:07:26 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn6.nearest_tss"

#' locusdef.rn6.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' Built on Sat Mar 18 13:07:26 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn6.nearest_gene"

#' locusdef.rn6.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' Built on Sat Mar 18 13:07:28 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn6.exon"

#' locusdef.rn6.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' Built on Sat Mar 18 13:07:30 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn6.intron"

#' locusdef.rn6.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 13:07:30 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn6.1kb"

#' locusdef.rn6.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 13:07:30 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn6.1kb_outside_upstream"

#' locusdef.rn6.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 13:07:31 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn6.1kb_outside"

#' locusdef.rn6.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 13:07:31 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn6.5kb"

#' locusdef.rn6.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 13:07:31 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn6.5kb_outside_upstream"

#' locusdef.rn6.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 13:07:32 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn6.5kb_outside"

#' locusdef.rn6.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' Built on Sat Mar 18 13:07:32 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn6.10kb"

#' locusdef.rn6.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Sat Mar 18 13:07:32 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn6.10kb_outside_upstream"

#' locusdef.rn6.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Sat Mar 18 13:07:33 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.0 and org.Rn.eg.db_3.4.0.
"locusdef.rn6.10kb_outside"

#' tss.dm3 TSS locations
#'
#' A \code{GRanges} with all the TSSs for dm3. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.4.0.
"tss.dm3"

#' locusdef.dm3.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:07:59 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.4.0.
"locusdef.dm3.nearest_tss"

#' locusdef.dm3.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:07:59 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.4.0.
"locusdef.dm3.nearest_gene"

#' locusdef.dm3.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:00 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.4.0.
"locusdef.dm3.exon"

#' locusdef.dm3.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:01 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.4.0.
"locusdef.dm3.intron"

#' locusdef.dm3.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:01 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.4.0.
"locusdef.dm3.1kb"

#' locusdef.dm3.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:02 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.4.0.
"locusdef.dm3.1kb_outside_upstream"

#' locusdef.dm3.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:02 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.4.0.
"locusdef.dm3.1kb_outside"

#' locusdef.dm3.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:02 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.4.0.
"locusdef.dm3.5kb"

#' locusdef.dm3.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:02 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.4.0.
"locusdef.dm3.5kb_outside_upstream"

#' locusdef.dm3.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:02 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.4.0.
"locusdef.dm3.5kb_outside"

#' locusdef.dm3.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:03 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.4.0.
"locusdef.dm3.10kb"

#' locusdef.dm3.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:03 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.4.0.
"locusdef.dm3.10kb_outside_upstream"

#' locusdef.dm3.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:03 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.4.0.
"locusdef.dm3.10kb_outside"

#' tss.dm6 TSS locations
#'
#' A \code{GRanges} with all the TSSs for dm6. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.3.0 and org.Dm.eg.db_3.4.0.
"tss.dm6"

#' locusdef.dm6.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:37 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.3.0 and org.Dm.eg.db_3.4.0.
"locusdef.dm6.nearest_tss"

#' locusdef.dm6.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:37 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.3.0 and org.Dm.eg.db_3.4.0.
"locusdef.dm6.nearest_gene"

#' locusdef.dm6.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:38 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.3.0 and org.Dm.eg.db_3.4.0.
"locusdef.dm6.exon"

#' locusdef.dm6.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:39 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.3.0 and org.Dm.eg.db_3.4.0.
"locusdef.dm6.intron"

#' locusdef.dm6.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:39 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.3.0 and org.Dm.eg.db_3.4.0.
"locusdef.dm6.1kb"

#' locusdef.dm6.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:39 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.3.0 and org.Dm.eg.db_3.4.0.
"locusdef.dm6.1kb_outside_upstream"

#' locusdef.dm6.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:39 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.3.0 and org.Dm.eg.db_3.4.0.
"locusdef.dm6.1kb_outside"

#' locusdef.dm6.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:40 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.3.0 and org.Dm.eg.db_3.4.0.
"locusdef.dm6.5kb"

#' locusdef.dm6.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:40 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.3.0 and org.Dm.eg.db_3.4.0.
"locusdef.dm6.5kb_outside_upstream"

#' locusdef.dm6.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:40 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.3.0 and org.Dm.eg.db_3.4.0.
"locusdef.dm6.5kb_outside"

#' locusdef.dm6.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:40 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.3.0 and org.Dm.eg.db_3.4.0.
"locusdef.dm6.10kb"

#' locusdef.dm6.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:40 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.3.0 and org.Dm.eg.db_3.4.0.
"locusdef.dm6.10kb_outside_upstream"

#' locusdef.dm6.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Sat Mar 18 13:08:40 2017.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.3.0 and org.Dm.eg.db_3.4.0.
"locusdef.dm6.10kb_outside"

#' geneset.GOBP.hsa genesets for Homo sapiens
#'
#' Gene Ontology Biological Process (GOBP) genesets for Homo sapiens. All genesets are required to >= 10 Entrez IDs.
#' Built on Sat Mar 18 13:14:28 2017.
#'
#' @format A \code{GeneSet} object with the following slots:
#' \describe{
#'     \item{type}{A \code{character} indicating the type of genesets, e.g. GOBP.}
#'     \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
#'     \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
#'     \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. GO IDs for GO terms), and whose elements are \code{character} vectors of Entrez Gene IDs.}
#'     \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
#'     \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
#' }
#' @source org.Hs.eg.db_3.4.0 and GO.db_3.4.0
"geneset.GOBP.hsa"

#' geneset.GOCC.hsa genesets for Homo sapiens
#'
#' Gene Ontology Cellular Component (GOCC) genesets for Homo sapiens. All genesets are required to >= 10 Entrez IDs.
#' Built on Sat Mar 18 13:14:28 2017.
#'
#' @format A \code{GeneSet} object with the following slots:
#' \describe{
#'     \item{type}{A \code{character} indicating the type of genesets, e.g. GOBP.}
#'     \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
#'     \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
#'     \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. GO IDs for GO terms), and whose elements are \code{character} vectors of Entrez Gene IDs.}
#'     \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
#'     \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
#' }
#' @source org.Hs.eg.db_3.4.0 and GO.db_3.4.0
"geneset.GOCC.hsa"

#' geneset.GOMF.hsa genesets for Homo sapiens
#'
#' Gene Ontology Molecular Function (GOMF) genesets for Homo sapiens. All genesets are required to >= 10 Entrez IDs.
#' Built on Sat Mar 18 13:14:28 2017.
#'
#' @format A \code{GeneSet} object with the following slots:
#' \describe{
#'     \item{type}{A \code{character} indicating the type of genesets, e.g. GOBP.}
#'     \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
#'     \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
#'     \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. GO IDs for GO terms), and whose elements are \code{character} vectors of Entrez Gene IDs.}
#'     \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
#'     \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
#' }
#' @source org.Hs.eg.db_3.4.0 and GO.db_3.4.0
"geneset.GOMF.hsa"

#' geneset.GOBP.mmu genesets for Mus musculus
#'
#' Gene Ontology Biological Process (GOBP) genesets for Mus musculus. All genesets are required to >= 10 Entrez IDs.
#' Built on Sat Mar 18 13:15:40 2017.
#'
#' @format A \code{GeneSet} object with the following slots:
#' \describe{
#'     \item{type}{A \code{character} indicating the type of genesets, e.g. GOBP.}
#'     \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
#'     \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
#'     \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. GO IDs for GO terms), and whose elements are \code{character} vectors of Entrez Gene IDs.}
#'     \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
#'     \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
#' }
#' @source org.Mm.eg.db_3.4.0 and GO.db_3.4.0
"geneset.GOBP.mmu"

#' geneset.GOCC.mmu genesets for Mus musculus
#'
#' Gene Ontology Cellular Component (GOCC) genesets for Mus musculus. All genesets are required to >= 10 Entrez IDs.
#' Built on Sat Mar 18 13:15:40 2017.
#'
#' @format A \code{GeneSet} object with the following slots:
#' \describe{
#'     \item{type}{A \code{character} indicating the type of genesets, e.g. GOBP.}
#'     \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
#'     \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
#'     \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. GO IDs for GO terms), and whose elements are \code{character} vectors of Entrez Gene IDs.}
#'     \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
#'     \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
#' }
#' @source org.Mm.eg.db_3.4.0 and GO.db_3.4.0
"geneset.GOCC.mmu"

#' geneset.GOMF.mmu genesets for Mus musculus
#'
#' Gene Ontology Molecular Function (GOMF) genesets for Mus musculus. All genesets are required to >= 10 Entrez IDs.
#' Built on Sat Mar 18 13:15:40 2017.
#'
#' @format A \code{GeneSet} object with the following slots:
#' \describe{
#'     \item{type}{A \code{character} indicating the type of genesets, e.g. GOBP.}
#'     \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
#'     \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
#'     \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. GO IDs for GO terms), and whose elements are \code{character} vectors of Entrez Gene IDs.}
#'     \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
#'     \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
#' }
#' @source org.Mm.eg.db_3.4.0 and GO.db_3.4.0
"geneset.GOMF.mmu"

#' geneset.GOBP.rno genesets for Rattus norvegicus
#'
#' Gene Ontology Biological Process (GOBP) genesets for Rattus norvegicus. All genesets are required to >= 10 Entrez IDs.
#' Built on Sat Mar 18 13:16:58 2017.
#'
#' @format A \code{GeneSet} object with the following slots:
#' \describe{
#'     \item{type}{A \code{character} indicating the type of genesets, e.g. GOBP.}
#'     \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
#'     \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
#'     \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. GO IDs for GO terms), and whose elements are \code{character} vectors of Entrez Gene IDs.}
#'     \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
#'     \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
#' }
#' @source org.Rn.eg.db_3.4.0 and GO.db_3.4.0
"geneset.GOBP.rno"

#' geneset.GOCC.rno genesets for Rattus norvegicus
#'
#' Gene Ontology Cellular Component (GOCC) genesets for Rattus norvegicus. All genesets are required to >= 10 Entrez IDs.
#' Built on Sat Mar 18 13:16:58 2017.
#'
#' @format A \code{GeneSet} object with the following slots:
#' \describe{
#'     \item{type}{A \code{character} indicating the type of genesets, e.g. GOBP.}
#'     \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
#'     \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
#'     \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. GO IDs for GO terms), and whose elements are \code{character} vectors of Entrez Gene IDs.}
#'     \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
#'     \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
#' }
#' @source org.Rn.eg.db_3.4.0 and GO.db_3.4.0
"geneset.GOCC.rno"

#' geneset.GOMF.rno genesets for Rattus norvegicus
#'
#' Gene Ontology Molecular Function (GOMF) genesets for Rattus norvegicus. All genesets are required to >= 10 Entrez IDs.
#' Built on Sat Mar 18 13:16:58 2017.
#'
#' @format A \code{GeneSet} object with the following slots:
#' \describe{
#'     \item{type}{A \code{character} indicating the type of genesets, e.g. GOBP.}
#'     \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
#'     \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
#'     \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. GO IDs for GO terms), and whose elements are \code{character} vectors of Entrez Gene IDs.}
#'     \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
#'     \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
#' }
#' @source org.Rn.eg.db_3.4.0 and GO.db_3.4.0
"geneset.GOMF.rno"

#' geneset.GOBP.dme genesets for Drosophila melanogaster
#'
#' Gene Ontology Biological Process (GOBP) genesets for Drosophila melanogaster. All genesets are required to >= 10 Entrez IDs.
#' Built on Sat Mar 18 13:17:29 2017.
#'
#' @format A \code{GeneSet} object with the following slots:
#' \describe{
#'     \item{type}{A \code{character} indicating the type of genesets, e.g. GOBP.}
#'     \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
#'     \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
#'     \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. GO IDs for GO terms), and whose elements are \code{character} vectors of Entrez Gene IDs.}
#'     \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
#'     \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
#' }
#' @source org.Dm.eg.db_3.4.0 and GO.db_3.4.0
"geneset.GOBP.dme"

#' geneset.GOCC.dme genesets for Drosophila melanogaster
#'
#' Gene Ontology Cellular Component (GOCC) genesets for Drosophila melanogaster. All genesets are required to >= 10 Entrez IDs.
#' Built on Sat Mar 18 13:17:29 2017.
#'
#' @format A \code{GeneSet} object with the following slots:
#' \describe{
#'     \item{type}{A \code{character} indicating the type of genesets, e.g. GOBP.}
#'     \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
#'     \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
#'     \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. GO IDs for GO terms), and whose elements are \code{character} vectors of Entrez Gene IDs.}
#'     \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
#'     \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
#' }
#' @source org.Dm.eg.db_3.4.0 and GO.db_3.4.0
"geneset.GOCC.dme"

#' geneset.GOMF.dme genesets for Drosophila melanogaster
#'
#' Gene Ontology Molecular Function (GOMF) genesets for Drosophila melanogaster. All genesets are required to >= 10 Entrez IDs.
#' Built on Sat Mar 18 13:17:29 2017.
#'
#' @format A \code{GeneSet} object with the following slots:
#' \describe{
#'     \item{type}{A \code{character} indicating the type of genesets, e.g. GOBP.}
#'     \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
#'     \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
#'     \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. GO IDs for GO terms), and whose elements are \code{character} vectors of Entrez Gene IDs.}
#'     \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
#'     \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
#' }
#' @source org.Dm.eg.db_3.4.0 and GO.db_3.4.0
"geneset.GOMF.dme"
#' geneset.reactome.hsa genesets for Homo sapiens
#'
#' Reactome genesets for Homo sapiens. All genesets are required to have >= 10 Entrez IDs.
#' Built on Sun Mar 19 16:11:23 2017.
#'
#' @format A \code{GeneSet} object with the following slots:
#' \describe{
#'     \item{type}{A \code{character} indicating the type of genesets, e.g. Reactome.}
#'     \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
#'     \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
#'     \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. R-HSA-109688), and whose elements are \code{character} vectors of Entrez Gene IDs.}
#'     \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
#'     \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
#' }
#' @source http://www.reactome.org/download/current/NCBI2Reactome_All_Levels.txt downloaded on 2017-03-19
"geneset.reactome.hsa"

#' geneset.reactome.mmu genesets for Mus musculus
#'
#' Reactome genesets for Mus musculus. All genesets are required to have >= 10 Entrez IDs.
#' Built on Sun Mar 19 16:11:28 2017.
#'
#' @format A \code{GeneSet} object with the following slots:
#' \describe{
#'     \item{type}{A \code{character} indicating the type of genesets, e.g. Reactome.}
#'     \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
#'     \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
#'     \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. R-HSA-109688), and whose elements are \code{character} vectors of Entrez Gene IDs.}
#'     \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
#'     \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
#' }
#' @source http://www.reactome.org/download/current/NCBI2Reactome_All_Levels.txt downloaded on 2017-03-19
"geneset.reactome.mmu"

#' geneset.reactome.rno genesets for Rattus norvegicus
#'
#' Reactome genesets for Rattus norvegicus. All genesets are required to have >= 10 Entrez IDs.
#' Built on Sun Mar 19 16:11:32 2017.
#'
#' @format A \code{GeneSet} object with the following slots:
#' \describe{
#'     \item{type}{A \code{character} indicating the type of genesets, e.g. Reactome.}
#'     \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
#'     \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
#'     \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. R-HSA-109688), and whose elements are \code{character} vectors of Entrez Gene IDs.}
#'     \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
#'     \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
#' }
#' @source http://www.reactome.org/download/current/NCBI2Reactome_All_Levels.txt downloaded on 2017-03-19
"geneset.reactome.rno"

#' geneset.reactome.dme genesets for Drosophila melanogaster
#'
#' Reactome genesets for Drosophila melanogaster. All genesets are required to have >= 10 Entrez IDs.
#' Built on Sun Mar 19 16:11:36 2017.
#'
#' @format A \code{GeneSet} object with the following slots:
#' \describe{
#'     \item{type}{A \code{character} indicating the type of genesets, e.g. Reactome.}
#'     \item{dburl}{A \code{character} of the URL of the database underlying the genesets.}
#'     \item{organism}{A \code{character} of the organism, e.g. Homo sapiens.}
#'     \item{set.gene}{An \code{environment} containing a \code{list} whose keys are database specific accessions (e.g. R-HSA-109688), and whose elements are \code{character} vectors of Entrez Gene IDs.}
#'     \item{all.genes}{A \code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \code{type}.}
#'     \item{set.name}{An \code{environment} containing a \code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}
#' }
#' @source http://www.reactome.org/download/current/NCBI2Reactome_All_Levels.txt downloaded on 2017-03-19
"geneset.reactome.dme"

