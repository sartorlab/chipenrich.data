#' tss.hg19 TSS locations
#'
#' A \code{GRanges} with all the TSSs for hg19. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.}
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism.
"tss.hg19"

#' locusdef.hg19.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg19.nearest_tss"

#' locusdef.hg19.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg19.nearest_gene"

#' locusdef.hg19.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg19.exon"

#' locusdef.hg19.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg19.intron"

#' locusdef.hg19.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg19.1kb"

#' locusdef.hg19.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg19.1kb_outside_upstream"

#' locusdef.hg19.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg19.1kb_outside"

#' locusdef.hg19.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg19.5kb"

#' locusdef.hg19.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg19.5kb_outside_upstream"

#' locusdef.hg19.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg19.5kb_outside"

#' locusdef.hg19.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg19.10kb"

#' locusdef.hg19.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg19.10kb_outside_upstream"

#' locusdef.hg19.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg19.10kb_outside"

#' tss.hg38 TSS locations
#'
#' A \code{GRanges} with all the TSSs for hg38. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.}
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism.
"tss.hg38"

#' locusdef.hg38.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg38.nearest_tss"

#' locusdef.hg38.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg38.nearest_gene"

#' locusdef.hg38.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg38.exon"

#' locusdef.hg38.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg38.intron"

#' locusdef.hg38.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg38.1kb"

#' locusdef.hg38.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg38.1kb_outside_upstream"

#' locusdef.hg38.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg38.1kb_outside"

#' locusdef.hg38.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg38.5kb"

#' locusdef.hg38.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg38.5kb_outside_upstream"

#' locusdef.hg38.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg38.5kb_outside"

#' locusdef.hg38.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg38.10kb"

#' locusdef.hg38.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg38.10kb_outside_upstream"

#' locusdef.hg38.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.hg38.10kb_outside"

#' tss.mm9 TSS locations
#'
#' A \code{GRanges} with all the TSSs for mm9. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.}
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism.
"tss.mm9"

#' locusdef.mm9.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm9.nearest_tss"

#' locusdef.mm9.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm9.nearest_gene"

#' locusdef.mm9.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm9.exon"

#' locusdef.mm9.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm9.intron"

#' locusdef.mm9.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm9.1kb"

#' locusdef.mm9.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm9.1kb_outside_upstream"

#' locusdef.mm9.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm9.1kb_outside"

#' locusdef.mm9.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm9.5kb"

#' locusdef.mm9.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm9.5kb_outside_upstream"

#' locusdef.mm9.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm9.5kb_outside"

#' locusdef.mm9.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm9.10kb"

#' locusdef.mm9.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm9.10kb_outside_upstream"

#' locusdef.mm9.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm9.10kb_outside"

#' tss.mm10 TSS locations
#'
#' A \code{GRanges} with all the TSSs for mm10. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.}
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism.
"tss.mm10"

#' locusdef.mm10.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm10.nearest_tss"

#' locusdef.mm10.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm10.nearest_gene"

#' locusdef.mm10.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm10.exon"

#' locusdef.mm10.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm10.intron"

#' locusdef.mm10.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm10.1kb"

#' locusdef.mm10.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm10.1kb_outside_upstream"

#' locusdef.mm10.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm10.1kb_outside"

#' locusdef.mm10.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm10.5kb"

#' locusdef.mm10.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm10.5kb_outside_upstream"

#' locusdef.mm10.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm10.5kb_outside"

#' locusdef.mm10.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm10.10kb"

#' locusdef.mm10.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm10.10kb_outside_upstream"

#' locusdef.mm10.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.mm10.10kb_outside"

#' tss.rn4 TSS locations
#'
#' A \code{GRanges} with all the TSSs for rn4. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.}
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism.
"tss.rn4"

#' locusdef.rn4.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn4.nearest_tss"

#' locusdef.rn4.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn4.nearest_gene"

#' locusdef.rn4.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn4.exon"

#' locusdef.rn4.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn4.intron"

#' locusdef.rn4.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn4.1kb"

#' locusdef.rn4.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn4.1kb_outside_upstream"

#' locusdef.rn4.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn4.1kb_outside"

#' locusdef.rn4.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn4.5kb"

#' locusdef.rn4.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn4.5kb_outside_upstream"

#' locusdef.rn4.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn4.5kb_outside"

#' locusdef.rn4.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn4.10kb"

#' locusdef.rn4.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn4.10kb_outside_upstream"

#' locusdef.rn4.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn4.10kb_outside"

#' tss.rn5 TSS locations
#'
#' A \code{GRanges} with all the TSSs for rn5. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.}
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism.
"tss.rn5"

#' locusdef.rn5.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn5.nearest_tss"

#' locusdef.rn5.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn5.nearest_gene"

#' locusdef.rn5.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn5.exon"

#' locusdef.rn5.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn5.intron"

#' locusdef.rn5.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn5.1kb"

#' locusdef.rn5.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn5.1kb_outside_upstream"

#' locusdef.rn5.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn5.1kb_outside"

#' locusdef.rn5.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn5.5kb"

#' locusdef.rn5.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn5.5kb_outside_upstream"

#' locusdef.rn5.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn5.5kb_outside"

#' locusdef.rn5.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn5.10kb"

#' locusdef.rn5.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn5.10kb_outside_upstream"

#' locusdef.rn5.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn5.10kb_outside"

#' tss.rn6 TSS locations
#'
#' A \code{GRanges} with all the TSSs for rn6. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.}
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism.
"tss.rn6"

#' locusdef.rn6.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn6.nearest_tss"

#' locusdef.rn6.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn6.nearest_gene"

#' locusdef.rn6.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn6.exon"

#' locusdef.rn6.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn6.intron"

#' locusdef.rn6.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn6.1kb"

#' locusdef.rn6.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn6.1kb_outside_upstream"

#' locusdef.rn6.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn6.1kb_outside"

#' locusdef.rn6.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn6.5kb"

#' locusdef.rn6.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn6.5kb_outside_upstream"

#' locusdef.rn6.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn6.5kb_outside"

#' locusdef.rn6.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn6.10kb"

#' locusdef.rn6.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn6.10kb_outside_upstream"

#' locusdef.rn6.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.rn6.10kb_outside"

#' tss.dm3 TSS locations
#'
#' A \code{GRanges} with all the TSSs for dm3. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.}
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism.
"tss.dm3"

#' locusdef.dm3.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm3.nearest_tss"

#' locusdef.dm3.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm3.nearest_gene"

#' locusdef.dm3.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm3.exon"

#' locusdef.dm3.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm3.intron"

#' locusdef.dm3.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm3.1kb"

#' locusdef.dm3.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm3.1kb_outside_upstream"

#' locusdef.dm3.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm3.1kb_outside"

#' locusdef.dm3.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm3.5kb"

#' locusdef.dm3.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm3.5kb_outside_upstream"

#' locusdef.dm3.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm3.5kb_outside"

#' locusdef.dm3.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm3.10kb"

#' locusdef.dm3.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm3.10kb_outside_upstream"

#' locusdef.dm3.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm3.10kb_outside"

#' tss.dm6 TSS locations
#'
#' A \code{GRanges} with all the TSSs for dm6. Primarily used in the \code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.}
#'
#' @format A \code{GRanges} object with the following \code{mcols}:
#' \describe{
#'     \item{gene_id}{The Entrez ID for the TSS}
#'     \item{symbol}{The gene symbol for the TSS}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism.
"tss.dm6"

#' locusdef.dm6.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm6.nearest_tss"

#' locusdef.dm6.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm6.nearest_gene"

#' locusdef.dm6.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm6.exon"

#' locusdef.dm6.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm6.intron"

#' locusdef.dm6.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm6.1kb"

#' locusdef.dm6.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm6.1kb_outside_upstream"

#' locusdef.dm6.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm6.1kb_outside"

#' locusdef.dm6.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm6.5kb"

#' locusdef.dm6.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm6.5kb_outside_upstream"

#' locusdef.dm6.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm6.5kb_outside"

#' locusdef.dm6.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm6.10kb"

#' locusdef.dm6.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm6.10kb_outside_upstream"

#' locusdef.dm6.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb_outside upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source Appropriate TxDb and orgDb packages for genome build and organism, and GENCODE Comprehensive gene annotations for reference chromosomes if available.
"locusdef.dm6.10kb_outside"

