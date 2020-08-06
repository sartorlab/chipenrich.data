#' locusdef.hg19.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' Built on Fri Apr 13 09:45:51 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.nearest_tss"

"locusdef.hg19.enhancer_plus5kb"
"locusdef.hg19.enhancer"
#' locusdef.hg19.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' Built on Fri Apr 13 09:45:51 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.nearest_gene"

#' locusdef.hg19.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' Built on Fri Apr 13 09:45:53 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.exon"

#' locusdef.hg19.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' Built on Fri Apr 13 09:45:55 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.intron"

#' locusdef.hg19.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:45:55 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.1kb"

#' locusdef.hg19.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:45:56 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.1kb_outside_upstream"

#' locusdef.hg19.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:45:56 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.1kb_outside"

#' locusdef.hg19.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:45:57 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.5kb"

#' locusdef.hg19.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:45:57 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.5kb_outside_upstream"

#' locusdef.hg19.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:45:58 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.5kb_outside"

#' locusdef.hg19.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:45:58 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.10kb"

#' locusdef.hg19.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:45:58 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.10kb_outside_upstream"

#' locusdef.hg19.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:45:59 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg19.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz
"locusdef.hg19.10kb_outside"

#' locusdef.hg38.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' Built on Fri Apr 13 09:49:34 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.nearest_tss"

#' locusdef.hg38.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' Built on Fri Apr 13 09:49:34 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.nearest_gene"

#' locusdef.hg38.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' Built on Fri Apr 13 09:49:37 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.exon"

#' locusdef.hg38.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' Built on Fri Apr 13 09:49:38 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.intron"

#' locusdef.hg38.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:49:39 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.1kb"

#' locusdef.hg38.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:49:39 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.1kb_outside_upstream"

#' locusdef.hg38.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:49:40 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.1kb_outside"

#' locusdef.hg38.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:49:41 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.5kb"

#' locusdef.hg38.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:49:41 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.5kb_outside_upstream"

#' locusdef.hg38.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:49:42 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.5kb_outside"

#' locusdef.hg38.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:49:42 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.10kb"

#' locusdef.hg38.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:49:43 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.10kb_outside_upstream"

#' locusdef.hg38.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:49:43 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, hg38.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Homo sapiens.}
#' }
#' @source R packages: TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0 and org.Hs.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz
"locusdef.hg38.10kb_outside"

#' locusdef.mm9.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' Built on Fri Apr 13 09:50:37 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.nearest_tss"

#' locusdef.mm9.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' Built on Fri Apr 13 09:50:37 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.nearest_gene"

#' locusdef.mm9.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' Built on Fri Apr 13 09:50:39 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.exon"

#' locusdef.mm9.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' Built on Fri Apr 13 09:50:41 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.intron"

#' locusdef.mm9.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:50:41 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.1kb"

#' locusdef.mm9.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:50:41 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.1kb_outside_upstream"

#' locusdef.mm9.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:50:42 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.1kb_outside"

#' locusdef.mm9.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:50:42 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.5kb"

#' locusdef.mm9.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:50:43 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.5kb_outside_upstream"

#' locusdef.mm9.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:50:43 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.5kb_outside"

#' locusdef.mm9.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:50:44 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.10kb"

#' locusdef.mm9.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:50:44 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.10kb_outside_upstream"

#' locusdef.mm9.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:50:44 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm9.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz
"locusdef.mm9.10kb_outside"

#' locusdef.mm10.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' Built on Fri Apr 13 09:51:33 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.nearest_tss"

#' locusdef.mm10.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' Built on Fri Apr 13 09:51:33 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.nearest_gene"

#' locusdef.mm10.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' Built on Fri Apr 13 09:51:35 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.exon"

#' locusdef.mm10.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' Built on Fri Apr 13 09:51:36 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.intron"

#' locusdef.mm10.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:51:37 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.1kb"

#' locusdef.mm10.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:51:37 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.1kb_outside_upstream"

#' locusdef.mm10.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:51:38 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.1kb_outside"

#' locusdef.mm10.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:51:38 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.5kb"

#' locusdef.mm10.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:51:39 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.5kb_outside_upstream"

#' locusdef.mm10.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:51:39 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.5kb_outside"

#' locusdef.mm10.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:51:40 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.10kb"

#' locusdef.mm10.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:51:40 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.10kb_outside_upstream"

#' locusdef.mm10.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:51:41 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, mm10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Mus musculus.}
#' }
#' @source R packages: TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0 and org.Mm.eg.db_3.5.0. GENCODE resources: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz and ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz
"locusdef.mm10.10kb_outside"

#' locusdef.rn4.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:51:58 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn4.nearest_tss"

#' locusdef.rn4.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:51:58 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn4.nearest_gene"

#' locusdef.rn4.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:51:59 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn4.exon"

#' locusdef.rn4.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:52:00 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn4.intron"

#' locusdef.rn4.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:52:01 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn4.1kb"

#' locusdef.rn4.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:52:01 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn4.1kb_outside_upstream"

#' locusdef.rn4.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:52:01 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn4.1kb_outside"

#' locusdef.rn4.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:52:01 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn4.5kb"

#' locusdef.rn4.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:52:02 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn4.5kb_outside_upstream"

#' locusdef.rn4.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:52:02 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn4.5kb_outside"

#' locusdef.rn4.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:52:02 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn4.10kb"

#' locusdef.rn4.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:52:02 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn4.10kb_outside_upstream"

#' locusdef.rn4.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the rn4 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:52:03 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn4.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn4.10kb_outside"

#' locusdef.rn5.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' Built on Fri Apr 13 09:52:16 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn5.nearest_tss"

#' locusdef.rn5.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' Built on Fri Apr 13 09:52:16 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn5.nearest_gene"

#' locusdef.rn5.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' Built on Fri Apr 13 09:52:17 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn5.exon"

#' locusdef.rn5.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' Built on Fri Apr 13 09:52:19 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn5.intron"

#' locusdef.rn5.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:52:19 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn5.1kb"

#' locusdef.rn5.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:52:19 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn5.1kb_outside_upstream"

#' locusdef.rn5.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:52:19 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn5.1kb_outside"

#' locusdef.rn5.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:52:20 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn5.5kb"

#' locusdef.rn5.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:52:20 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn5.5kb_outside_upstream"

#' locusdef.rn5.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:52:20 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn5.5kb_outside"

#' locusdef.rn5.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:52:20 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn5.10kb"

#' locusdef.rn5.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:52:21 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn5.10kb_outside_upstream"

#' locusdef.rn5.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:52:21 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn5.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn5.refGene_3.4.2 and org.Rn.eg.db_3.5.0.
"locusdef.rn5.10kb_outside"

#' locusdef.rn6.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' Built on Fri Apr 13 09:52:34 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.1 and org.Rn.eg.db_3.5.0.
"locusdef.rn6.nearest_tss"

#' locusdef.rn6.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' Built on Fri Apr 13 09:52:35 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.1 and org.Rn.eg.db_3.5.0.
"locusdef.rn6.nearest_gene"

#' locusdef.rn6.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' Built on Fri Apr 13 09:52:36 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.1 and org.Rn.eg.db_3.5.0.
"locusdef.rn6.exon"

#' locusdef.rn6.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' Built on Fri Apr 13 09:52:37 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.1 and org.Rn.eg.db_3.5.0.
"locusdef.rn6.intron"

#' locusdef.rn6.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:52:38 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.1 and org.Rn.eg.db_3.5.0.
"locusdef.rn6.1kb"

#' locusdef.rn6.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:52:38 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.1 and org.Rn.eg.db_3.5.0.
"locusdef.rn6.1kb_outside_upstream"

#' locusdef.rn6.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:52:38 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.1 and org.Rn.eg.db_3.5.0.
"locusdef.rn6.1kb_outside"

#' locusdef.rn6.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:52:38 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.1 and org.Rn.eg.db_3.5.0.
"locusdef.rn6.5kb"

#' locusdef.rn6.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:52:39 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.1 and org.Rn.eg.db_3.5.0.
"locusdef.rn6.5kb_outside_upstream"

#' locusdef.rn6.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:52:39 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.1 and org.Rn.eg.db_3.5.0.
"locusdef.rn6.5kb_outside"

#' locusdef.rn6.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:52:39 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.1 and org.Rn.eg.db_3.5.0.
"locusdef.rn6.10kb"

#' locusdef.rn6.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:52:39 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.1 and org.Rn.eg.db_3.5.0.
"locusdef.rn6.10kb_outside_upstream"

#' locusdef.rn6.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:52:40 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, rn6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Rattus norvegicus.}
#' }
#' @source R packages: TxDb.Rnorvegicus.UCSC.rn6.refGene_3.4.1 and org.Rn.eg.db_3.5.0.
"locusdef.rn6.10kb_outside"

#' locusdef.dm3.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:03 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.5.0.
"locusdef.dm3.nearest_tss"

#' locusdef.dm3.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:03 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.5.0.
"locusdef.dm3.nearest_gene"

#' locusdef.dm3.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:04 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.5.0.
"locusdef.dm3.exon"

#' locusdef.dm3.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:04 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.5.0.
"locusdef.dm3.intron"

#' locusdef.dm3.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:05 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.5.0.
"locusdef.dm3.1kb"

#' locusdef.dm3.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:05 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.5.0.
"locusdef.dm3.1kb_outside_upstream"

#' locusdef.dm3.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:05 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.5.0.
"locusdef.dm3.1kb_outside"

#' locusdef.dm3.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:05 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.5.0.
"locusdef.dm3.5kb"

#' locusdef.dm3.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:05 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.5.0.
"locusdef.dm3.5kb_outside_upstream"

#' locusdef.dm3.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:05 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.5.0.
"locusdef.dm3.5kb_outside"

#' locusdef.dm3.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:05 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.5.0.
"locusdef.dm3.10kb"

#' locusdef.dm3.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:05 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.5.0.
"locusdef.dm3.10kb_outside_upstream"

#' locusdef.dm3.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm3 genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:06 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm3.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 and org.Dm.eg.db_3.5.0.
"locusdef.dm3.10kb_outside"

#' locusdef.dm6.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:33 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.4.1 and org.Dm.eg.db_3.5.0.
"locusdef.dm6.nearest_tss"

#' locusdef.dm6.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:33 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.4.1 and org.Dm.eg.db_3.5.0.
"locusdef.dm6.nearest_gene"

#' locusdef.dm6.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:33 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.4.1 and org.Dm.eg.db_3.5.0.
"locusdef.dm6.exon"

#' locusdef.dm6.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:34 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.4.1 and org.Dm.eg.db_3.5.0.
"locusdef.dm6.intron"

#' locusdef.dm6.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:34 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.4.1 and org.Dm.eg.db_3.5.0.
"locusdef.dm6.1kb"

#' locusdef.dm6.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:34 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.4.1 and org.Dm.eg.db_3.5.0.
"locusdef.dm6.1kb_outside_upstream"

#' locusdef.dm6.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:35 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.4.1 and org.Dm.eg.db_3.5.0.
"locusdef.dm6.1kb_outside"

#' locusdef.dm6.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:35 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.4.1 and org.Dm.eg.db_3.5.0.
"locusdef.dm6.5kb"

#' locusdef.dm6.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:35 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.4.1 and org.Dm.eg.db_3.5.0.
"locusdef.dm6.5kb_outside_upstream"

#' locusdef.dm6.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:35 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.4.1 and org.Dm.eg.db_3.5.0.
"locusdef.dm6.5kb_outside"

#' locusdef.dm6.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:35 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.4.1 and org.Dm.eg.db_3.5.0.
"locusdef.dm6.10kb"

#' locusdef.dm6.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:35 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.4.1 and org.Dm.eg.db_3.5.0.
"locusdef.dm6.10kb_outside_upstream"

#' locusdef.dm6.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' For the dm6 genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.
#'
#' Built on Fri Apr 13 09:53:35 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, dm6.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Drosophila melanogaster.}
#' }
#' @source R packages: TxDb.Dmelanogaster.UCSC.dm6.ensGene_3.4.1 and org.Dm.eg.db_3.5.0.
"locusdef.dm6.10kb_outside"

#' locusdef.danRer10.nearest_tss locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs.
#'
#' Built on Fri Apr 13 09:54:22 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, danRer10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Danio rerio.}
#' }
#' @source R packages: TxDb.Drerio.UCSC.danRer10.refGene_3.4.2 and org.Dr.eg.db_3.5.0.
"locusdef.danRer10.nearest_tss"

#' locusdef.danRer10.nearest_gene locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs.
#'
#' Built on Fri Apr 13 09:54:23 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, danRer10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Danio rerio.}
#' }
#' @source R packages: TxDb.Drerio.UCSC.danRer10.refGene_3.4.2 and org.Dr.eg.db_3.5.0.
"locusdef.danRer10.nearest_gene"

#' locusdef.danRer10.exon locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.
#'
#' Built on Fri Apr 13 09:54:24 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, danRer10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Danio rerio.}
#' }
#' @source R packages: TxDb.Drerio.UCSC.danRer10.refGene_3.4.2 and org.Dr.eg.db_3.5.0.
"locusdef.danRer10.exon"

#' locusdef.danRer10.intron locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.
#'
#' Built on Fri Apr 13 09:54:25 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, danRer10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Danio rerio.}
#' }
#' @source R packages: TxDb.Drerio.UCSC.danRer10.refGene_3.4.2 and org.Dr.eg.db_3.5.0.
"locusdef.danRer10.intron"

#' locusdef.danRer10.1kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 1kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:54:25 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, danRer10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Danio rerio.}
#' }
#' @source R packages: TxDb.Drerio.UCSC.danRer10.refGene_3.4.2 and org.Dr.eg.db_3.5.0.
"locusdef.danRer10.1kb"

#' locusdef.danRer10.1kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:54:25 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, danRer10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Danio rerio.}
#' }
#' @source R packages: TxDb.Drerio.UCSC.danRer10.refGene_3.4.2 and org.Dr.eg.db_3.5.0.
"locusdef.danRer10.1kb_outside_upstream"

#' locusdef.danRer10.1kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 1kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:54:26 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, danRer10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Danio rerio.}
#' }
#' @source R packages: TxDb.Drerio.UCSC.danRer10.refGene_3.4.2 and org.Dr.eg.db_3.5.0.
"locusdef.danRer10.1kb_outside"

#' locusdef.danRer10.5kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 5kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:54:26 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, danRer10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Danio rerio.}
#' }
#' @source R packages: TxDb.Drerio.UCSC.danRer10.refGene_3.4.2 and org.Dr.eg.db_3.5.0.
"locusdef.danRer10.5kb"

#' locusdef.danRer10.5kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:54:26 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, danRer10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Danio rerio.}
#' }
#' @source R packages: TxDb.Drerio.UCSC.danRer10.refGene_3.4.2 and org.Dr.eg.db_3.5.0.
"locusdef.danRer10.5kb_outside_upstream"

#' locusdef.danRer10.5kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 5kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:54:27 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, danRer10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Danio rerio.}
#' }
#' @source R packages: TxDb.Drerio.UCSC.danRer10.refGene_3.4.2 and org.Dr.eg.db_3.5.0.
"locusdef.danRer10.5kb_outside"

#' locusdef.danRer10.10kb locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined within 10kb upstream and downstream of the TSS.
#'
#' Built on Fri Apr 13 09:54:27 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, danRer10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Danio rerio.}
#' }
#' @source R packages: TxDb.Drerio.UCSC.danRer10.refGene_3.4.2 and org.Dr.eg.db_3.5.0.
"locusdef.danRer10.10kb"

#' locusdef.danRer10.10kb_outside_upstream locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.
#'
#' Built on Fri Apr 13 09:54:27 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, danRer10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Danio rerio.}
#' }
#' @source R packages: TxDb.Drerio.UCSC.danRer10.refGene_3.4.2 and org.Dr.eg.db_3.5.0.
"locusdef.danRer10.10kb_outside_upstream"

#' locusdef.danRer10.10kb_outside locus definition
#'
#' A \code{LocusDefinition} where a gene locus is defined as the region beyond 10kb upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.
#'
#' Built on Fri Apr 13 09:54:27 2018.
#'
#' @format A \code{LocusDefinition} object with the following slots:
#' \describe{
#'     \item{granges}{A \code{GRanges} of the locus definitions with \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}}
#'     \item{dframe}{A \code{data.frame} of the locus definitions with columns for \code{chr}, \code{start}, \code{end}, \code{gene_id}, and \code{symbol}}
#'     \item{genome.build}{A \code{character} indicating the genome build. In this case, danRer10.}
#'     \item{organism}{A \code{character} indicating the organism name. In this case, Danio rerio.}
#' }
#' @source R packages: TxDb.Drerio.UCSC.danRer10.refGene_3.4.2 and org.Dr.eg.db_3.5.0.
"locusdef.danRer10.10kb_outside"

