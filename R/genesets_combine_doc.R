#' geneset.biocarta_pathway.hsa genesets for BioCarta
#'
#' BioCarta (biocarta_pathway) genesets. All genesets are required to have >= 10 Entrez IDs.
#' Built on Mon Oct 16 18:46:04 2017.
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
#' @source https://cgap.nci.nih.gov/Pathways/BioCarta_Pathways
"geneset.biocarta_pathway.hsa"

#' geneset.ctd.hsa genesets for Comparative Toxicogenomics Database
#'
#' Comparative Toxicogenomics Database (ctd) genesets. All genesets are required to have >= 10 Entrez IDs.
#' Built on Mon Oct 16 18:46:11 2017.
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
#' @source http://ctdbase.org
"geneset.ctd.hsa"

#' geneset.drug_bank.hsa genesets for DrugBank
#'
#' DrugBank (drug_bank) genesets. All genesets are required to have >= 10 Entrez IDs.
#' Built on Mon Oct 16 18:46:13 2017.
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
#' @source https://www.drugbank.ca
"geneset.drug_bank.hsa"

#' geneset.hallmark.hsa genesets for Hallmark (MSigDB)
#'
#' Hallmark (MSigDB) (hallmark) genesets. All genesets are required to have >= 10 Entrez IDs.
#' Built on Mon Oct 16 18:46:15 2017.
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
#' @source http://software.broadinstitute.org/gsea/msigdb/collections.jsp#H
"geneset.hallmark.hsa"

#' geneset.immunologic.hsa genesets for Immunologic Signatures (MSigDB)
#'
#' Immunologic Signatures (MSigDB) (immunologic) genesets. All genesets are required to have >= 10 Entrez IDs.
#' Built on Mon Oct 16 18:46:45 2017.
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
#' @source http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C7
"geneset.immunologic.hsa"

#' geneset.kegg_pathway.hsa genesets for KEGG Pathways
#'
#' KEGG Pathways (kegg_pathway) genesets. All genesets are required to have >= 10 Entrez IDs.
#' Built on Mon Oct 16 18:46:53 2017.
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
#' @source http://kegg.jp
"geneset.kegg_pathway.hsa"

#' geneset.microrna.hsa genesets for MicroRNA Targets (MSigDB)
#'
#' MicroRNA Targets (MSigDB) (microrna) genesets. All genesets are required to have >= 10 Entrez IDs.
#' Built on Mon Oct 16 18:46:56 2017.
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
#' @source http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C3
"geneset.microrna.hsa"

#' geneset.oncogenic.hsa genesets for Oncogenic Signatures (MSigDB)
#'
#' Oncogenic Signatures (MSigDB) (oncogenic) genesets. All genesets are required to have >= 10 Entrez IDs.
#' Built on Mon Oct 16 18:47:23 2017.
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
#' @source http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C6
"geneset.oncogenic.hsa"

#' geneset.pfam.hsa genesets for Pfam
#'
#' Pfam (pfam) genesets. All genesets are required to have >= 10 Entrez IDs.
#' Built on Mon Oct 16 18:47:28 2017.
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
#' @source http://pfam.xfam.org
"geneset.pfam.hsa"

#' geneset.transcription_factors.hsa genesets for Transcription Factor Targets (MSigDB)
#'
#' Transcription Factor Targets (MSigDB) (transcription_factors) genesets. All genesets are required to have >= 10 Entrez IDs.
#' Built on Mon Oct 16 18:47:33 2017.
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
#' @source http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C3
"geneset.transcription_factors.hsa"

#' geneset.protein_interaction_biogrid.hsa genesets for BioGRID Protein Interactions
#'
#' BioGRID Protein Interactions (protein_interaction_biogrid) genesets. All genesets are required to have >= 10 Entrez IDs.
#' Built on Tue Oct 24 16:05:53 2017.
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
#' @source https://thebiogrid.org
"geneset.protein_interaction_biogrid.hsa"

