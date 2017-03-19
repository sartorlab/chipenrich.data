get_package_version = function(package_name) {
	packages = data.frame(installed.packages()[,c('Package','Version')], row.names = NULL, stringsAsFactors=F)
	version_str = paste(subset(packages, Package == package_name), collapse='_')
	return(version_str)
}

# The precede() and follow() functions are quiet weird so it is recommended that
# you look at a simple example and try to figure out what's going on.
#
# In particular, precede() and follow() return indices in subject for which the
# element of x is the direct precedent or follower, respectively, of. It would
# appear that precede() and follow() get it backwards.
#
# a = GenomicRanges::GRanges(
#     seqnames = c('chr1','chr1','chr1'),
#     ranges = IRanges::IRanges(start = c(1,3,5), end = c(1,3,5)),
#     strand = '*')
#
# b = GenomicRanges::GRanges(
#     seqnames = c('chr1','chr1','chr1'),
#     ranges = IRanges::IRanges(start = c(2,4,6), end = c(2,4,6)),
#     strand = '*')
#
# precede(x = a, subject = a)
# follow(a, a)
#
# precede(a, b)
# follow(a, b)
#
# precede(b, a)
# follow(b, a)

################################################################################
### Workflow
###                                txdb                    orgdb
###                                  |                       |
###                                  -------------------------
###                                  |           |           |
###                  build_transcript_gr()   build_exons_introns_gr()
###                                  |           |           |
###                        gr_transcripts     gr_exons     gr_introns
###                                  |           |           |
###                                  -------------------------
###                                              |
###                  gencode + ens2eg -> filter_gencode_types()
###                                              |
###                  eg2symbol -> filter_readthrough_transcripts()
###                                              |
###                                  -------------------------
###                                  |           |           |
###                        gr_transcripts     gr_exons   gr_introns
###                                  |           |           |
###                 add_ntss_precursors()        |           |
###                                  |      reduce_gr()  reduce_gr()
###                add_ngene_precursors()        |           |
###                                  |           |           |
###       add_nkb_precursors(width = 1000)       |           |
###       add_nkb_precursors(width = 5000)       |           |
###      add_nkb_precursors(width = 10000)       |           |
###                                  |           |           |
###               build_ldef_nearest_tss()       |           |
###                                  |           |           |
###              build_ldef_nearest_gene()       |           |
###                                  |           |           |
###          build_ldef_nkb(1000,'inside')       |           |
###        build_ldef_nkb(1000,'upstream')       |           |
###         build_ldef_nkb(1000,'outside')       |           |
###                                  |           |           |
###          build_ldef_nkb(5000,'inside')       |           |
###        build_ldef_nkb(5000,'upstream')       |           |
###         build_ldef_nkb(5000,'outside')       |           |
###                                  |           |           |
###         build_ldef_nkb(10000,'inside')       |           |
###       build_ldef_nkb(10000,'upstream')       |           |
###        build_ldef_nkb(10000,'outside')       |           |
###                                  |           |           |
###                                  -------------------------
###                                              |
###                               convert to LocusDefinition class
###                                              |
###                                       output to /data
###
################################################################################

################################################################################
### Build base GRanges
################################################################################

###
### Transcripts
###
#' Function to build initial transcript GRanges for locus definitions
#'
#' A function to build initial transcript \code{GRanges} that are used to build
#' nearest_tss, nearest_gene, and nkb-related locus definitions. If the \code{TxDb}
#' object uses ENSEMBL (rn4 or dm3) or FLYBASE IDs (dm6), the \code{ensembl2eg}
#' parameter must be used to convert them to Entrez IDs to get gene symbols.
#'
#' Transcripts without Entrez IDs are removed because genesets are in terms of Entrez IDs.
#'
#' @param txdb A TxDb object as from the numerous TxDb.* packages in Bioconductor.
#' @param eg2symbol A data.frame created from the \code{egSYMBOL} object in the \code{orgDb} packages. Required.
#' @param ensembl2eg A data.frame created from the \code{egENSEMBL2EG} object in the \code{orgDb} packages. Use \code{NULL} if gene IDs are not ENSEMBL or FLYBASE.
#'
#' @return A sorted \code{GRanges} with Entrez ID and symbol metadata that is unique up to location and Entrez ID.
#'

build_transcript_gr = function(txdb, eg2symbol, ensembl2eg) {
    if(is.null(ensembl2eg)) {
        # Get transcripts with Entrez IDs and add gene symbols
        gr = GenomicFeatures::transcripts(txdb, columns=c('gene_id'))
        gr$gene_id = as.integer(gr$gene_id)
        gr$symbol = eg2symbol[match(gr$gene_id, eg2symbol$gene_id), 'symbol']
    } else {
        # Get transcripts with ENSEMBL IDs
        gr = GenomicFeatures::transcripts(txdb, columns=c('gene_id'))
        gr$gene_id = as.character(gr$gene_id)

        if(unique(genome(gr)) == 'dm6') {
            # Then these are not actually ENSEMBL IDs, but FLYBASE IDs ... UGH
            # And we need to strip the .1 from the names ... UGH!
            flybase_geneid_list = strsplit(gr$gene_id, '[.]')
            gr$gene_id = sapply(flybase_geneid_list, function(l){l[1]})
        }

        # Rename gene_id mcol to ensembl_id
        colnames(mcols(gr)) = 'ensembl_id'
        # Get Entrez IDs from ENSEMBL IDs
        gr$gene_id = ensembl2eg[match(gr$ensembl_id, ensembl2eg$ensembl_id), 'gene_id']
        # Get symbols from Entrez IDs
        gr$symbol = eg2symbol[match(gr$gene_id, eg2symbol$gene_id), 'symbol']
    }

    # Keep only gene_id and symbol columns
    GenomicRanges::mcols(gr) = GenomicRanges::mcols(gr)[,c('gene_id','symbol')]

    # Force the range to have a gene_id
    # Can't do chipenrich without it
    gr = gr[!is.na(gr$gene_id)]

    # Enforce uniqueness on location + gene_id
    gr = unique(gr)
    gr = sort(gr)

    return(gr)
}

###
### Exons and introns
###
#' Function to build GRanges exon and intron locus definitions
#'
#' A function to build initial exon/intron \code{GRanges} that are used to build
#' nearest_tss, nearest_gene, and nkb-related locus definitions. If the \code{TxDb}
#' object uses ENSEMBL (rn4 or dm3) or FLYBASE IDs (dm6), the \code{ensembl2eg}
#' parameter must be used to convert them to Entrez IDs to get gene symbols.
#'
#' Transcripts without Entrez IDs are removed because genesets are in terms of Entrez IDs.
#'
#' @param txdb A TxDb object as from the numerous TxDb.* packages in Bioconductor.
#' @param eg2symbol A data.frame created from the \code{egSYMBOL} object in the \code{orgDb} packages.
#' @param ensembl2eg An optional data.frame created from the \code{egENSEMBL2EG} object in the \code{orgDb} packages.
#' @param type One of 'exons' or 'introns'.
#'
#' @return A sorted \code{GRanges} with Entrez ID and symbol metadata that is unique up to location and Entrez ID.
#'

build_exons_introns_gr = function(txdb, eg2symbol, ensembl2eg, type = c('exons','introns')) {

    # NOTE: For rn5, rn6 there is a warning:
    # Warning message:
    # In .set_group_names(grl, use.names, txdb, by) :
    #   some group names are NAs or duplicated
    # This occurs because the exon_name column has NAs this is not a problem.
    if(type == 'exons') {
        grl = GenomicFeatures::exonsBy(txdb, by= 'tx', use.names = TRUE)
    } else {
        grl = GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)
    }

    # This is a trick to quickly construct txname vectors for the GRangesList
    # which has the txname in the name() slot of the list, but not in the GRanges
    # contained in the list
    txname_rle = S4Vectors::Rle(names(grl), S4Vectors::elementNROWS(grl))
    txname_vec = as.character(txname_rle)

    gr = unlist(grl, use.names = FALSE)
    mcols(gr)$tx_name = txname_vec

    if(is.null(ensembl2eg)) {
        # Get mapping from UCSC tx_name to Entrez Gene ID
        # UCSC TXID and TXNAME to GENEID mapping (for introns and exons)
        id_maps = AnnotationDbi::select(txdb, keys = names(grl), columns = c('TXID','GENEID'), keytype = 'TXNAME')
        # Add Entrez Gene ID from UCSC tx_name
        mcols(gr)$gene_id = id_maps[match(mcols(gr)$tx_name, id_maps$TXNAME), 'GENEID']
        # Add symbol from Entrez Gene ID
        mcols(gr)$symbol = eg2symbol[match(mcols(gr)$gene_id, eg2symbol$gene_id), 'symbol']
    } else {
        # Get the ENSEMBL Transcript ID to ENSEMBL Gene ID mapping
        # When TxDb uses ENSEMBL, the result of exonsBy() or intronsByTranscript()
        # gives ENSEMBL Transcript IDs in the tx_name column of the GRanges
        id_maps = AnnotationDbi::select(txdb, keys = names(grl), columns = c('TXID','GENEID'), keytype = 'TXNAME')
        mcols(gr)$ensembl_id = id_maps[match(mcols(gr)$tx_name, id_maps$TXNAME), 'GENEID']

        if(unique(genome(gr)) == 'dm6') {
            # Then these are not actually ENSEMBL IDs, but FLYBASE IDs ... UGH
            # And we need to strip the .1 from the names ... UGH!
            flybase_geneid_list = strsplit(gr$ensembl_id, '[.]')
            gr$ensembl_id = sapply(flybase_geneid_list, function(l){l[1]})
        }

        # Add Entrez Gene ID from ENSEMBL Gene ID
        mcols(gr)$gene_id = ensembl2eg[match(gr$ensembl_id, ensembl2eg$ensembl_id), 'gene_id']
        # Add symbol from Entrez Gene ID
        mcols(gr)$symbol = eg2symbol[match(mcols(gr)$gene_id, eg2symbol$gene_id), 'symbol']
    }

    # Keep only gene_id and symbol in mcols
    mcols(gr) = mcols(gr)[,c('gene_id','symbol')]

    # Force the gene_id to integer from CharacterList
    # NOTE: Worth a check that mean(sapply(gr$gene_id, length)) = 1
    # so we're sure a transcript doesn't have multiple gene_ids.
    gr$gene_id = as.integer(gr$gene_id)

    # Force the range to have a gene_id
    # Can't do chipenrich without it
    gr = gr[!is.na(gr$gene_id)]

    # Enforce uniqueness on location + gene_id
    gr = unique(gr)
    gr = sort(gr)

    return(gr)
}

################################################################################
### Filter Entrez Gene IDs that correspond to certain GENCODE types
################################################################################

#' Function to filter certain GENCODE biotypes
#'
#' After creating initial \code{GRanges} with \code{build_transcript_gr()} or
#' \code{build_exons_introns_gr()}, the Entrez IDs belonging to certain GENCODE
#' biotypes (https://www.gencodegenes.org/gencode_biotypes.html) are removed.
#'
#' In particular, biotypes filters include: '3prime_overlapping_ncRNA', 'scRNA',
#' 'snoRNA', 'snRNA', 'scaRNA', anti-sense transcripts, and those beginning with 'LOC'.
#'
#' If a genome does not have annotations from GENCODE, i.e. rat and fly, nothing
#' is done, but this can be changed at a later point.
#'
#' @param gr A \code{GRanges} resulting from \code{build_transcript_gr()} or \code{build_exons_introns_gr()}.
#' @param gencode A \code{GRanges} of the GENCODE Comprehensive gene annotation on reference chromosomes.
#' @param ens2eg A \code{data.frame} mapping Entrez Ids to ENSEMBL IDs (identifiers in \code{gencode}).
#'
#' @return A \code{GRanges} with listed biotypes filtered out.
#'

filter_gencode_types = function(gr, gencode, ens2eg) {
    # If we have gencode data, use it for filtering
    if(!is.null(gencode) && !is.null(ens2eg)) {
        # gene_type codes from GENCODE to filter out
        filter_gencode_types = c('3prime_overlapping_ncRNA', 'scRNA', 'snoRNA', 'snRNA', 'scaRNA')

        filter_types_ensembl_ids = unique(unlist(subset(gencode, gene_type %in% filter_gencode_types)$Parent, use.names = FALSE))

        # Filter gene_name with -AS for antisense
        filter_AS_ensembl_ids = unique(unlist(gencode[grepl('-AS', gencode$gene_name), 'Parent'], use.names = FALSE))

        # Combine
        filter_ensembl_ids = c(filter_types_ensembl_ids, filter_AS_ensembl_ids)

        # Get the Entrez IDs to use in gr, gr_exons, and gr_introns
        filter_entrez_ids = unique(unlist(subset(ens2eg, ens_id %in% filter_ensembl_ids)$gene_id, use.names = FALSE))

        # Filter from gr
        gr = gr[!(gr$gene_id %in% filter_entrez_ids)]

        # Filter from gr those with symbols beginning with LOC
        gr = gr[!grepl('^LOC', gr$symbol)]
    } else {
        # Not sure what to do yet
    }

    return(gr)
}

################################################################################
### Filter readthrough symbols, e.g. APITD1-CORT
################################################################################

#' Function to filter features from readthroughs
#'
#' After filtering with \code{filter_gencode_types()}, \code{GRanges} are filtered
#' if they belong to a readthrough gene, i.e. APITD1-CORT.
#'
#' We also look for gene symbols having a '-' where both strings on either side
#' are valid gene symbols. We use the \code{egGENENAME} object in the \code{orgDb}
#' packages, to look for genes with 'through' in the name. This deals with the
#' few cases where aliases are used and thus missed with the '-' strategy.
#'
#' @param gr A \code{GRanges} after going through \code{filter_gencode_types()}.
#' @param egGENENAME A data.frame created from the \code{egGENENAME} object in the \code{orgDb} packages.
#'
#' @return A \code{GRanges} bereft of features in readthrough genes.
#'

filter_readthrough_transcripts = function(gr, egGENENAME) {
    ### Based on 'through' in gene name
        # Get a data.frame of Entrez IDs mapped to gene names
        mapped_genes = AnnotationDbi::mappedkeys(egGENENAME)
        gene_names = as.data.frame(egGENENAME[mapped_genes])

        # Gene Entrez IDs whose gene names have 'through' in the name
        through_idx = grep('through', gene_names$gene_name)
        through_egid = unique(gene_names[through_idx, 'gene_id'])

    ### Based on - alone
        dash_readthroughs_idx = grep('-', gr$symbol)
        dash_readthroughs = mcols(gr[dash_readthroughs_idx])[, c('gene_id','symbol')]

        # Take this as the universe of possible
        symbol_universe = grep('-', gr$symbol, value=T, invert=T)

        # Split the possible_readthroughs by the -
        dash_readthroughs_split = strsplit(dash_readthroughs$symbol, '-')

        # Figure out which are valid
        is_dash_readthrough = sapply(dash_readthroughs_split, function(drs){
            return(all(drs %in% symbol_universe))
        })

        dash_readthroughs = dash_readthroughs[is_dash_readthrough,]
        dash_readthroughs_egid = unique(dash_readthroughs$gene_id)

    ### Combine the Entrez IDs from both methods
        combined_readthroughs_egid = unique(c(through_egid, dash_readthroughs_egid))

    ### Get rid of the things in gr with Entrez IDs in combined_readthroughs_egid
    gr = gr[!(gr$gene_id %in% combined_readthroughs_egid)]

    return(gr)
}

################################################################################
### Tweaks to precede() and follow() to fix when either function returns NA.
### Turn NA into the end of the chromosome (follow) or the start (precede).
### Use the GenomeInfoDb::seqinfo() of the subject GRanges to determine length of chromosomes.
### NOTE: GenomeInfoDb::seqinfo(x) should be the same as GenomeInfoDb::seqinfo(subject).
################################################################################

#' Helper function to find preceding \code{GRanges}
#'
#' A beefed up wrapper of the \code{GenomicRanges::precede} function, that uses
#' the \code{seqinfo} of \code{subject} to add chromosome starts, so no NAs result.
#'
#' An addition fix is required when any feature occurs exactly at the start of
#' a chromosome. In such cases, +5bp is added to the appropriate ranges in \code{x}.
#'
#' NOTE: Internally, we actually use \code{GenomicRanges::follow} to get the preceding
#' item. The documentation is a bit confusing, but this does the correct thing.
#'
#' @param x A \code{GRanges} object having point ranges.
#' @param subject A \code{GRanges} object having point ranges.
#'
#' @return A \code{GRanges} with same length as \code{x} that has the range in
#' \code{subject} that precedes the corresponding range for the row in \code{x}.
#' This function is used to add a \code{GRanges} to \code{x}.
#'

subject_preceding_x = function(x, subject) {
    # Create GRanges of starts of chromosomes
    gr_chr_start = GenomicRanges::GRanges(
        seqnames = as.character(GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(subject))),
        ranges = IRanges::IRanges(start = 1, end = 1),
        strand = '*',
        seqinfo = GenomeInfoDb::seqinfo(subject))

    # Find overlaps between x and gr_chr_start and shift the start of those by 5bp
    # This prevents genes that start or end at the chromosome boundary from causing problems
    x_overlaps_start = findOverlaps(x, gr_chr_start, ignore.strand = TRUE)
    if(length(x_overlaps_start) > 0) {
        start(x)[queryHits(x_overlaps_start)] = start(x)[queryHits(x_overlaps_start)] + 5
        end(x)[queryHits(x_overlaps_start)] = end(x)[queryHits(x_overlaps_start)] + 5
    }

    # Add start ranges to subject so there are no NAs in follow()
    subject = sort(c(subject, gr_chr_start))

    # Determine indices of ranges in subject that precede items in x
    precede_idx = follow(x, subject, ignore.strand = TRUE)

    # Create initial draft of GRanges to return
    gr = granges(subject[precede_idx])

    return(gr)
}

#' Helper function to find following \code{GRanges}
#'
#' A beefed up wrapper of the \code{GenomicRanges::follow} function, that uses
#' the \code{seqinfo} of \code{subject} to add chromosome ends, so no NAs result.
#'
#' An addition fix is required when any feature occurs exactly at the end of
#' a chromosome. In such cases, -5bp is subjtracted to the appropriate ranges in \code{x}.
#'
#' NOTE: Internally, we actually use \code{GenomicRanges::follow} to get the preceding
#' item. The documentation is a bit confusing, but this does the correct thing.
#'
#' @param x A \code{GRanges} object having point ranges.
#' @param subject A \code{GRanges} object having point ranges.
#'
#' @return A \code{GRanges} with same length as \code{x} that has the range in
#' \code{subject} that follows the corresponding range for the row in \code{x}.
#' This function is used to add a \code{GRanges} to \code{x}.
#'

subject_following_x = function(x, subject) {
    # Create GRanges of starts of chromosomes
    lengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(subject))
    gr_chr_end = GenomicRanges::GRanges(
        seqnames = as.character(GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(subject))),
        ranges = IRanges::IRanges(start = lengths, end = lengths),
        strand = '*',
        seqinfo = GenomeInfoDb::seqinfo(subject))

    # Find overlaps between x and gr_chr_end and shift the end of those by 5bp
    # This prevents genes that start or end at the chromosome boundary from causing problems
    x_overlaps_end = findOverlaps(x, gr_chr_end, ignore.strand = TRUE)
    if(length(x_overlaps_end) > 0) {
        start(x)[queryHits(x_overlaps_end)] = start(x)[queryHits(x_overlaps_end)] - 5
        end(x)[queryHits(x_overlaps_end)] = end(x)[queryHits(x_overlaps_end)] - 5
    }

    # Add end ranges to subject so there are no NAs in precede()
    subject = sort(c(subject, gr_chr_end))

    # Determine indices of ranges in subject that precede items in x
    follow_idx = precede(x, subject, ignore.strand = TRUE)

    # Create initial draft of GRanges to return
    gr = granges(subject[follow_idx])

    return(gr)
}

################################################################################
### Everything you need for TSS and nkb definitions right here!
################################################################################
#' Function to add required elements for nearest_tss locus definition
#'
#' The approach taken is to represent all TSSs as point \code{GRanges}, find the
#' preceding and following TSSs for each (using chromosome ends as needed), create
#' \code{IRanges} that span the space between TSSs, and use \code{mid} to set the
#' boundaries for the nearest_tss definition.
#'
#' @param gr A \code{GRanges} object after calling \code{build_transcript_gr()} and the filter functions.
#'
#' @return The \code{gr} \code{GRanges} with \code{mcols} for building the nearest_tss locus definition.
#'

add_ntss_precursors = function(gr) {
    pos_idx = which(strand(gr) == '+')
    neg_idx = which(strand(gr) == '-')

    # Add an integer column for the TSS
    mcols(gr)$tss = 0
    mcols(gr[pos_idx])$tss = start(gr)[pos_idx]
    mcols(gr[neg_idx])$tss = end(gr)[neg_idx]

    # Add a GRanges column for the TSS
    # NOTE: Use this for the tss.[genome_build].RData
    mcols(gr)$gr_tss = GenomicRanges::GRanges(
        seqnames = GenomeInfoDb::seqnames(gr),
        ranges = IRanges::IRanges(start = gr$tss, end = gr$tss),
        strand = strand(gr),
        seqinfo = GenomeInfoDb::seqinfo(gr))

    # Add GRanges columns for the following and preceding TSSs
    gr$tss_follow_tss = subject_following_x(gr$gr_tss, gr$gr_tss)
    gr$tss_precede_tss = subject_preceding_x(gr$gr_tss, gr$gr_tss)

    ### Determine boundaries of nearest_tss definitions. Essentially, midpoints btw TSSs

    # Add IRanges columns for starts and ends of locus def for transcripts
    # which will be used to take the midpoint
    gr$ntss_start_range = IRanges::IRanges(start = start(gr$tss_precede_tss), end = start(gr$gr_tss))
    gr$ntss_end_range = IRanges::IRanges(start = start(gr$gr_tss), end = start(gr$tss_follow_tss))

    # Add integer columns for start and end of locus def for transcript
    gr$ntss_start = mid(gr$ntss_start_range) + 1
    gr$ntss_end = mid(gr$ntss_end_range)

    # Add GRanges column of the final nearest_tss locus def
    gr$nearest_tss = GenomicRanges::GRanges(
        seqnames = GenomeInfoDb::seqnames(gr),
        ranges = IRanges::IRanges(start = gr$ntss_start, end = gr$ntss_end),
        strand = '*')

    return(gr)
}

################################################################################
### Add gr_left and gr_right mcols() to the input GRanges
################################################################################
#' Function to add required elements for nearest_gene locus definition
#'
#' The approach is to treat the left and right ends of the transcripts as we treated
#' the TSSs in \code{add_ntss_precursors()}. Create \code{GRanges} representing
#' the left and right ends of each row, find the preceding and following features
#' for each among the combined left and right ends, create \code{IRanges} that
#' extend from the preceding to the following item of each end of the feature,
#' use the midpoint to define the boundaries of the nearest_gene boundaries.
#'
#' NOTE: The nearest_gene locus definition is preferential to local information,
#' so genes contained in other genes remain whole and their containing genes are split.
#'
#' @param gr A \code{GRanges} object after calling \code{build_transcript_gr()}, the filter functions, and \code{add_ntss_precursors()}.
#'
#' @return The \code{gr} \code{GRanges} with additional \code{mcols} for building the nearest_gene locus definition.
#'

add_ngene_precursors = function(gr) {
    # Add GRanges for the left and right end of the feature
    mcols(gr)$gr_left = GenomicRanges::GRanges(
        seqnames = GenomeInfoDb::seqnames(gr),
        ranges = IRanges::IRanges(start = start(gr), end = start(gr)),
        strand = strand(gr),
        seqinfo = GenomeInfoDb::seqinfo(gr))

    mcols(gr)$gr_right = GenomicRanges::GRanges(
        seqnames = GenomeInfoDb::seqnames(gr),
        ranges = IRanges::IRanges(start = end(gr), end = end(gr)),
        strand = strand(gr),
        seqinfo = GenomeInfoDb::seqinfo(gr))

    # Combine the left and rights so it will find whichever is closest
    gr_combined = sort(c(gr$gr_right, gr$gr_left))

    # Add GRanges for precede and follow for left and right ends
    gr$precede_right = subject_preceding_x(gr$gr_right, gr_combined)
    gr$follow_right = subject_following_x(gr$gr_right, gr_combined)

    gr$precede_left = subject_preceding_x(gr$gr_left, gr_combined)
    gr$follow_left = subject_following_x(gr$gr_left, gr_combined)

    # Add IRanges with which to take midpoints for the nearest_gene definition
    ngene_left_start_range = IRanges::IRanges(start = start(gr$precede_left), end = start(gr))
    ngene_left_end_range = IRanges::IRanges(start = start(gr), end = start(gr$follow_left))

    ngene_right_start_range = IRanges::IRanges(start = start(gr$precede_right), end = end(gr))
    ngene_right_end_range = IRanges::IRanges(start = end(gr), end = start(gr$follow_right))

    # Add integer columns representing the start and end of the left and right sides
    # The nearest_gene boundaries are defined as the midpoints between the left
    # end of the range and the closest preceding TSS or TES and the right end of the
    # range and the closest following TSS or TES.
    gr$ngene_left_start = mid(ngene_left_start_range) + 1
    gr$ngene_left_end = mid(ngene_left_end_range)

    gr$ngene_right_start = mid(ngene_right_start_range) + 1
    gr$ngene_right_end = mid(ngene_right_end_range)

    return(gr)
}

################################################################################
### Create nkb related objects
################################################################################

#' Function to add required elements for nkb-related locus definitions
#'
#' The approach is to \code{GenomicRanges::flank} each TSS by \code{width} and then
#' restrict each flanked range by the preceding or following TSS.
#'
#' Also built in this function are regions outside of \code{width} and up to the
#' next TSS on just the upstream side, and on either side.
#'
#' @param gr A \code{GRanges} object after calling \code{build_transcript_gr()}, the filter functions, \code{add_ntss_precursors()}, and \code{add_ngene_precursors()}.
#' @param width An integer indicating the number of base-pairs for the locus definition. One of 1000, 5000, or 10000.
#'
#' @return The \code{gr} \code{GRanges} with additional \code{mcols} building within \code{width}kb, outside \code{width}kb and upstream, and outside \code{width}kb locus definitions.
#'

add_nkb_precursors = function(gr, width) {
    ### Create shortcodes for nkbs
    if(width == 1000) {
        short_code = '1kb'
    } else if (width == 5000) {
        short_code = '5kb'
    } else if (width == 10000) {
        short_code = '10kb'
    } else {
        stop('width should be one of 1000, 5000, 10000')
    }

    ### nkb base
        gr_nkb = flank(
            x = gr$gr_tss,
            width = width,
            both = TRUE,
            use.names = FALSE,
            ignore.strand = TRUE)
        GenomeInfoDb::seqinfo(gr_nkb) = GenomeInfoDb::seqinfo(gr)

    ### Inside nkb of TSS
        ### Create nkb definition by determining if the nkb start or end exceeds
        ### ntss_start and ntss_end, respectively.
        nkb_start_less_ntss_start_idx = which(start(gr_nkb) < gr$ntss_start)
        nkb_end_greater_ntss_end_idx = which(end(gr_nkb) > gr$ntss_end)

        # Fix the nkb ranges that go too far and add gene_id and symbol columns
        gr_nkb_fix = gr_nkb
        start(gr_nkb_fix)[nkb_start_less_ntss_start_idx] = gr$ntss_start[nkb_start_less_ntss_start_idx]
        end(gr_nkb_fix)[nkb_end_greater_ntss_end_idx] = gr$ntss_end[nkb_end_greater_ntss_end_idx]
        gr_nkb_fix$gene_id = gr$gene_id
        gr_nkb_fix$symbol = gr$symbol

    ### Outside nkb of TSS
        # Use restrict to get left and right side of 1kb def up to nearest_tss
        gr_outside_nkb_left = restrict(gr$nearest_tss, start = start(gr$nearest_tss), end = start(gr_nkb_fix))
        gr_outside_nkb_right = restrict(gr$nearest_tss, start = end(gr_nkb_fix), end = end(gr$nearest_tss))

        # Determine upstream ranges based on the strand
        gr_nkb_upstream_pos = gr_outside_nkb_left[strand(gr) == '+']
        mcols(gr_nkb_upstream_pos) = mcols(gr)[strand(gr) == '+', c('gene_id','symbol')]
        strand(gr_nkb_upstream_pos) = '+'

        gr_nkb_upstream_neg = gr_outside_nkb_right[strand(gr) == '-']
        mcols(gr_nkb_upstream_neg) = mcols(gr)[strand(gr) == '-', c('gene_id','symbol')]
        strand(gr_nkb_upstream_neg) = '-'

        # Determine downstream ranges based on the strand
        gr_nkb_downstream_pos = gr_outside_nkb_right[strand(gr) == '+']
        mcols(gr_nkb_downstream_pos) = mcols(gr)[strand(gr) == '+', c('gene_id','symbol')]
        strand(gr_nkb_downstream_pos) = '+'

        gr_nkb_downstream_neg = gr_outside_nkb_left[strand(gr) == '-']
        mcols(gr_nkb_downstream_neg) = mcols(gr)[strand(gr) == '-', c('gene_id','symbol')]
        strand(gr_nkb_downstream_neg) = '-'

        # Recombine the upstream and downstream strands
        gr_nkb_upstream = c(gr_nkb_upstream_pos, gr_nkb_upstream_neg)
        GenomeInfoDb::seqinfo(gr_nkb_upstream) = GenomeInfoDb::seqinfo(gr)

        gr_nkb_downstream = c(gr_nkb_downstream_pos, gr_nkb_downstream_neg)
        GenomeInfoDb::seqinfo(gr_nkb_downstream) = GenomeInfoDb::seqinfo(gr)

    ### Add gr_nkb_fix, gr_nkb_upstream, and gr_nkb_downstream to mcols(gr)
        gr$gr_nkb = trim(gr_nkb_fix)
        gr$gr_nkb_upstream = trim(gr_nkb_upstream)
        gr$gr_nkb_downstream = trim(gr_nkb_downstream)

        # Rename the columns
        colnames(mcols(gr))[which(colnames(mcols(gr)) == 'gr_nkb')] = sprintf('gr_%s', short_code)
        colnames(mcols(gr))[which(colnames(mcols(gr)) == 'gr_nkb_upstream')] = sprintf('gr_%s_upstream', short_code)
        colnames(mcols(gr))[which(colnames(mcols(gr)) == 'gr_nkb_downstream')] = sprintf('gr_%s_downstream', short_code)

    return(gr)
}

################################################################################
### A function to combine adjacent ranges belonging to the same gene_id so that
### locus definitions have the minimum number of lines possible.
################################################################################

#' Helper function to reduce overlapping regions with the same Entrez ID
#'
#' Used in all \code{build_ldef_*()} functions.
#'
#' @param gr A \code{GRanges} locus definition just prior to building \code{LocusDefinition} class.
#' @param eg2symbol A data.frame created from the \code{egSYMBOL} object in the \code{orgDb} packages.
#'
#' @return A \code{GRanges} where overlapping regions with the same Entrez ID are
#' combined, and the \code{gene_id} and \code{symbol} are preserved. This is the finishing
#' step prior to creating a \code{LocusDefinition} object.
#'

reduce_gr = function(gr, eg2symbol) {
    # Splitting on gene_id into GRangesList
    tmp_grl = splitAsList(gr, gr$gene_id)

    # Reducing within gene_id
    tmp_grl = reduce(tmp_grl)

    # Get the correct vector of the gene_id from the GRangesList
    gene_id_rle = S4Vectors::Rle(names(tmp_grl), S4Vectors::elementNROWS(tmp_grl))
    gene_id_vec = as.integer(gene_id_rle)

    # Reconstructing GRanges
    gr = unlist(tmp_grl, use.names = FALSE)

    # Appending gene_id and symbol since reduce() destroys all mcols()
    mcols(gr)$gene_id = gene_id_vec
    mcols(gr)$symbol = eg2symbol[match(mcols(gr)$gene_id, eg2symbol$gene_id), 'symbol']

    gr = sort(gr)

    return(gr)
}

################################################################################
################################################################################
################################################################################

################################################################################
### Build nearest_tss definition
################################################################################

#' Function to create the nearest_tss \code{GRanges}
#'
#' This function is called after \code{add_ntss_precursors()} adds all the precursors.
#' NOTE: At this step, regions with the same Entrez ID that overlap are reduced with
#' \code{reduce_gr()}.
#'
#' @param gr A \code{GRanges} locus definition just prior to building \code{LocusDefinition} class.
#' @param eg2symbol A data.frame created from the \code{egSYMBOL} object in the \code{orgDb} packages.
#'
#' @return A \code{GRanges} to be used in the \code{granges} slot of the nearest_tss \code{LocusDefinition} object.
#'

build_ldef_nearest_tss = function(gr, eg2symbol) {
    # Extract ntss
    ldef = gr$nearest_tss

    # Add gene_id and symbols
    ldef$gene_id = gr$gene_id
    ldef$symbol = gr$symbol

    # Inherit GenomeInfoDb::seqinfo()
    GenomeInfoDb::seqinfo(ldef) = GenomeInfoDb::seqinfo(gr)

    # Remove strand information
    strand(ldef) = '*'

    # Sort and unique
    ldef = sort(ldef)
    ldef = unique(ldef)

    # Reduce
    ldef = reduce_gr(ldef, eg2symbol)

    return(ldef)
}

################################################################################
### Build nearest_gene definition
################################################################################

#' Function to create the nearest_gene \code{GRanges}
#'
#' This function is called after \code{add_ngene_precursors()} adds all the precursors.
#' NOTE: At this step, regions with the same Entrez ID that overlap are reduced with
#' \code{reduce_gr()}.
#'
#' @param gr A \code{GRanges} locus definition just prior to building \code{LocusDefinition} class.
#' @param eg2symbol A data.frame created from the \code{egSYMBOL} object in the \code{orgDb} packages.
#'
#' @return A \code{GRanges} to be used in the \code{granges} slot of the nearest_gene \code{LocusDefinition} object.
#'

build_ldef_nearest_gene = function(gr, eg2symbol) {
    ### Create the nearest_gene definition by concatenating the GRanges consisting of:
    ### [gr$precede_left, start(gr)] and [start(gr), gr$follow_left]
    ### [gr$precede_right, end(gr)] and [end(gr), gr$follow_right]
    ### In other words, we are expanding around each transcript's left and right ends,
    ### ranges up until the preceding and following. There will initially be some
    ### redundancy that will be fixed with reduce() in a GRangesList on the gene_id.

    # Build the nearest_gene definition
    ldef = c(
        GenomicRanges::GRanges(
            seqnames = GenomeInfoDb::seqnames(gr),
            ranges = IRanges::IRanges(start = gr$ngene_left_start, end = gr$ngene_left_end),
            strand = '*',
            gene_id = gr$gene_id,
            symbol = gr$symbol),
        GenomicRanges::GRanges(
            seqnames = GenomeInfoDb::seqnames(gr),
            ranges = IRanges::IRanges(start = gr$ngene_right_start, end = gr$ngene_right_end),
            strand = '*',
            gene_id = gr$gene_id,
            symbol = gr$symbol)
    )

    # Enforce uniqueness of ranges/gene_id
    ldef = sort(ldef)
    ldef = unique(ldef)

    # Reduce
    ldef = reduce_gr(ldef, eg2symbol)

    return(ldef)
}

################################################################################
### Build nkb related definition
################################################################################

#' Function to create the nkb \code{GRanges}
#'
#' This function is called after \code{add_nkb_precursors()} adds all the precursors.
#' NOTE: At this step, regions with the same Entrez ID that overlap are reduced with
#' \code{reduce_gr()}.
#'
#' @param gr A \code{GRanges} locus definition just prior to building \code{LocusDefinition} class.
#' @param width One of 1000, 5000, or 10000 indicating how wide regions flanking TSSs should be.
#' @param context One of 'inside', 'upstream', or 'outside', indicating whether the regions should be
#' within, outside and upstream, or just outside of \code{width}.
#' @param eg2symbol A data.frame created from the \code{egSYMBOL} object in the \code{orgDb} packages.
#'
#' @return A \code{GRanges} to be used in the \code{granges} slot of the appropriate nkb \code{LocusDefinition} object.
#'

build_ldef_nkb = function(gr, width, context, eg2symbol) {

    # Convert to *kb
    if(width == 1000) {
        short_width = '1kb'
    } else if (width == 5000) {
        short_width = '5kb'
    } else if (width == 10000) {
        short_width = '10kb'
    }

    # Build the expected column names for extraction from mcols()
    if(context == 'inside') {
        short_codes = sprintf('gr_%s', short_width)
    } else if (context == 'upstream') {
        short_codes = sprintf('gr_%s_upstream', short_width)
    } else if (context == 'outside') {
        short_codes = c(
            sprintf('gr_%s_upstream', short_width),
            sprintf('gr_%s_downstream', short_width))
    }

    # Determine which
    mcols_idx =  colnames(mcols(gr)) %in% short_codes

    ldef = mcols(gr)[, mcols_idx]

    # If outside nkb, then concatenate the two resulting columns
    # otherwise, ldef is already just a GRanges
    if(context == 'outside') {
        ldef = c(ldef[,1], ldef[,2])
    }

    # Force unstranded
    strand(ldef) = '*'

    # Only take regions with width > 1
    ldef = ldef[which(width(ldef) > 1)]

    # Enforce uniqueness of ranges/gene_id
    ldef = sort(ldef)
    ldef = unique(ldef)

    # Reduce
    ldef = reduce_gr(ldef, eg2symbol)

    return(ldef)
}

################################################################################
### Build BED output
################################################################################

#' Function to write a \code{GRanges} locus definition as a BED file
#'
#' This function writes BED files for use in the Genome Browser to check locus definitions.
#'
#' @param ldef_gr A \code{GRanges} representing the final locus definition.
#' @param ldef_path A full path, including file name, indicating where to write BED file.
#'
#' @return \code{NULL}
#'

ldef_gr_to_bed_file = function(ldef_gr, ldef_path) {
    df = data.frame(ldef_gr, stringsAsFactors=F)
    df$strand = '.'
    df$score = 1000
    df$name = paste(df$gene_id, df$symbol, sep=':')

    df = df[,c('seqnames','start','end','name','score','strand')]
    write.table(df, file = ldef_path, sep='\t', quote = F, col.names=F, row.names=F)

    return(NULL)
}

# old_ldef_to_bed_file = function(old_ldef, ldef_name) {
#     df = data.frame(old_ldef@granges, stringsAsFactors=F)
#     df$strand = '.'
#     df$score = 1000
#
#     df = df[,c('seqnames','start','end','names','score','strand')]
#     write.table(df, file = sprintf('data/%s.bed', ldef_name), sep='\t', quote = F, col.names=F, row.names=F)
# }

#' Function to write a \code{GRanges} for the TSS RData
#'
#' This function creates the TSS.[genome_build] objects to be stored in \code{data/}
#'
#' @param gr A \code{GRanges} after calling \code{add_ntss_precursors}.
#' @param genome The genome build, e.g. hg19.
#'
#' @return \code{NULL}
#'

to_tss_rdata = function(gr, genome, txdb_version, orgdb_version, gencode_version, mapping_version) {
    tss = gr$gr_tss
    mcols(tss)$gene_id = mcols(gr)$gene_id
    mcols(tss)$symbol = mcols(gr)$symbol

    tss = unique(tss)
    tss = sort(tss)

    # The gr_tss object in mcols(gr) is exactly what's needed
    assign(sprintf('tss.%s', genome), tss)
    save(list = sprintf('tss.%s', genome), file = sprintf('data/tss.%s.RData', genome), compress = 'xz')

    ### Document the tss object
    # The actual documentation to write, as a roxygen2 block
    # See: http://r-pkgs.had.co.nz/data.html section on Documenting datasets

	# Additional information about GENCODE
	if(genome %in% c('hg19', 'hg38', 'mm9', 'mm10')) {
		source_desc = sprintf("#' @source R packages: %s and %s. GENCODE resources: %s and %s", txdb_version, orgdb_version, gencode_version, mapping_version)
	} else {
		source_desc = sprintf("#' @source R packages: %s and %s.", txdb_version, orgdb_version)
	}

    doc = c(
        sprintf("#' tss.%s TSS locations", genome),
        "#'",
        sprintf("#' A \\code{GRanges} with all the TSSs for %s. Primarily used in the \\code{assign_peaks()} function to report distance of a peak to the nearest TSS. Also used to build the QC plot with distribution of peaks to TSSs.", genome),
        "#'",
        "#' @format A \\code{GRanges} object with the following \\code{mcols}:",
        "#' \\describe{",
        "#'     \\item{gene_id}{The Entrez ID for the TSS}",
        "#'     \\item{symbol}{The gene symbol for the TSS}",
        "#' }",
        source_desc,
        sprintf('"tss.%s"', genome),
        ''
    )

    write(doc, file = 'R/data.R', append = TRUE)
}

#' Function to write a \code{GRanges} locus definition as a BED file
#'
#' This function writes \code{RData} \code{LocusDefinition} objects to \code{data/}.
#' Codes are of the form locusdef.[genome].[ldef_name].
#'
#' @param ldef_gr A \code{GRanges} representing the final locus definition.
#' @param genome The genome build, e.g. hg19.
#' @param organism Full name of organism, e.g. Homo sapiens.
#' @param ldef_name One of \code{c('nearest_tss', 'nearest_gene', 'exon', 'intron', '1kb', '1kb_outside_upstream', '1kb_outside', '5kb', '5kb_outside_upstream', '5kb_outside', '10kb', '10kb_outside_upstream','10kb_outside')}.
#'
#' @return \code{NULL}
#'

ldef_gr_to_LocusDefinition = function(ldef_gr, genome, organism, ldef_name, txdb_version, orgdb_version, gencode_version ,mapping_version) {

    # Convert to data.frame
    df = data.frame(ldef_gr, stringsAsFactors=FALSE)
    df = df[,c('seqnames','start','end','gene_id','symbol')]
    colnames(df) = c('chr','start','end','gene_id','symbol')

    # Build the LocusDefinition
    object = new("LocusDefinition")
    object@granges = ldef_gr
    object@dframe = df
    object@genome.build = genome
    object@organism = organism

    # Assign the locus definition to the proper name, and save it to data/
    assign(sprintf('locusdef.%s.%s', genome, ldef_name), object)
    save(list = sprintf('locusdef.%s.%s', genome, ldef_name), file = sprintf('data/locusdef.%s.%s.RData', genome, ldef_name), compress = 'xz')

    ### Document it and add the documentation to R/data.R

	# Detailed description based on specific locus definition
    if(ldef_name == 'nearest_tss') {
        detailed_desc = "#' A \\code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs."
    } else if (ldef_name == 'nearest_gene') {
        detailed_desc = "#' A \\code{LocusDefinition} where a gene locus is defined as the region spanning the midpoints between adjacent TSSs and TESs."
    } else if (ldef_name == '1kb' || ldef_name == '5kb' || ldef_name == '10kb') {
        detailed_desc = sprintf("#' A \\code{LocusDefinition} where a gene locus is defined within %s upstream and downstream of the TSS.", ldef_name)
    } else if (ldef_name == '1kb_outside_upstream' || ldef_name == '5kb_outside_upstream' || ldef_name == '10kb_outside_upstream') {
        nkb = unlist(strsplit(ldef_name, '_'))[1]
        detailed_desc = sprintf("#' A \\code{LocusDefinition} where a gene locus is defined as the region beyond %s upstream of the TSS and bounded by the midpoint between the TSS and the next upstream TSS.", nkb)
    } else if (ldef_name == '1kb_outside' || ldef_name == '5kb_outside' || ldef_name == '10kb_outside') {
        nkb = unlist(strsplit(ldef_name, '_'))[1]
        detailed_desc = sprintf("#' A \\code{LocusDefinition} where a gene locus is defined as the region beyond %s upstream and downstream of the TSS and bounded by the midpoints between the TSS and the next upstream and downstream TSSs.", nkb)
    } else if (ldef_name == 'exon') {
        detailed_desc = sprintf("#' A \\code{LocusDefinition} where a gene locus is defined as the exons belonging to genes.", ldef_name)
    } else if (ldef_name == 'intron') {
        detailed_desc = sprintf("#' A \\code{LocusDefinition} where a gene locus is defined as the introns belonging to genes.", ldef_name)
    }

    # Additional information based on the genome
    if(genome == 'rn4' || genome == 'dm3') {
        genome_desc = c(
            "#'",
            sprintf("#' For the %s genome, original gene IDs are from ENSEMBL and so an additional step of converting to Entrez IDs is done.", genome),
            "#'")
    } else if (genome == 'dm6') {
        genome_desc = c(
            "#'",
            sprintf("#' For the %s genome, original gene IDs are from FLYBASE and so an additional step of converting to Entrez IDs is done.", genome),
            "#'")
    } else {
        genome_desc = "#'"
    }

	# Additional information about GENCODE
	if(genome %in% c('hg19', 'hg38', 'mm9', 'mm10')) {
		source_desc = sprintf("#' @source R packages: %s and %s. GENCODE resources: %s and %s", txdb_version, orgdb_version, gencode_version, mapping_version)
	} else {
		source_desc = sprintf("#' @source R packages: %s and %s.", txdb_version, orgdb_version)
	}

    # The actual documentation to write, as a roxygen2 block
    # See: http://r-pkgs.had.co.nz/data.html section on Documenting datasets
    doc = c(
        sprintf("#' locusdef.%s.%s locus definition", genome, ldef_name),
        "#'",
        detailed_desc,
        genome_desc,
        sprintf("#' Built on %s.", date()),
        "#'",
        "#' @format A \\code{LocusDefinition} object with the following slots:",
        "#' \\describe{",
        "#'     \\item{granges}{A \\code{GRanges} of the locus definitions with \\code{mcols} for Entrez Gene ID \\code{gene_id} and gene symbol \\code{symbol}}",
        "#'     \\item{dframe}{A \\code{data.frame} of the locus definitions with columns for \\code{chr}, \\code{start}, \\code{end}, \\code{gene_id}, and \\code{symbol}}",
        sprintf("#'     \\item{genome.build}{A \\code{character} indicating the genome build. In this case, %s.}", genome),
        sprintf("#'     \\item{organism}{A \\code{character} indicating the organism name. In this case, %s.}", organism),
        "#' }",
        source_desc,
        sprintf('"locusdef.%s.%s"', genome, ldef_name),
        ''
    )

    write(doc, file = 'R/data.R', append = TRUE)

	return(NULL)
}

################################################################################
### The function to rule them all
################################################################################
#' Function to create all locus definitions for a given organism
#'
#' This function puts together all the other helper functions in the package to
#' build and write all locus definitions to \code{data/} as \code{RData}.
#'
#' @param genome One of 'hg19','hg38','mm9','mm10','rn4','rn5','rn6','dm3','dm6'
#'
#' @return \code{NULL}
#'

build_locus_definitions = function(genome) {

    ############################################################################
    ### Establish all required variables based on the genome
    ############################################################################

    if(genome == 'hg19') {
        # Gives Entrez IDs
        require(TxDb.Hsapiens.UCSC.hg19.knownGene)
        require(org.Hs.eg.db)
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
        egSYMBOL = org.Hs.egSYMBOL
        egENSEMBL2EG = NULL
        egGENENAME = org.Hs.egGENENAME
        gencode_url = 'data-raw/gencode.v25lift37.annotation.gff3.gz'
        mapping_url = 'data-raw/gencode.v25lift37.metadata.EntrezGene.gz'
        organism = 'Homo sapiens'
        txdb_version = get_package_version('TxDb.Hsapiens.UCSC.hg19.knownGene')
        orgdb_version = get_package_version('org.Hs.eg.db')
		gencode_version = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz'
		mapping_version = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz'
    } else if (genome == 'hg38') {
        # Gives Entrez IDs
        require(TxDb.Hsapiens.UCSC.hg38.knownGene)
        require(org.Hs.eg.db)
        txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
        egSYMBOL = org.Hs.egSYMBOL
        egENSEMBL2EG = NULL
        egGENENAME = org.Hs.egGENENAME
        gencode_url = 'data-raw/gencode.v25.annotation.gff3.gz'
        mapping_url = 'data-raw/gencode.v25.metadata.EntrezGene.gz'
        organism = 'Homo sapiens'
        txdb_version = get_package_version('TxDb.Hsapiens.UCSC.hg38.knownGene')
        orgdb_version = get_package_version('org.Hs.eg.db')
		gencode_version = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz'
		mapping_version = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.metadata.EntrezGene.gz'
    } else if (genome == 'mm9') {
        # Gives Entrez IDs
        require(TxDb.Mmusculus.UCSC.mm9.knownGene)
        require(org.Mm.eg.db)
        txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
        egSYMBOL = org.Mm.egSYMBOL
        egENSEMBL2EG = NULL
        egGENENAME = org.Mm.egGENENAME
        gencode_url = 'data-raw/gencode.vM9.annotation.gff3.gz'
        mapping_url = 'data-raw/gencode.vM9.metadata.EntrezGene.gz'
        organism = 'Mus musculus'
        txdb_version = get_package_version('TxDb.Mmusculus.UCSC.mm9.knownGene')
        orgdb_version = get_package_version('org.Mm.eg.db')
		gencode_version = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gff3.gz'
		mapping_version = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.metadata.EntrezGene.gz'
    } else if (genome == 'mm10') {
        # Gives Entrez IDs
        require(TxDb.Mmusculus.UCSC.mm10.knownGene)
        require(org.Mm.eg.db)
        txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
        egSYMBOL = org.Mm.egSYMBOL
        egENSEMBL2EG = NULL
        egGENENAME = org.Mm.egGENENAME
        gencode_url = 'data-raw/gencode.vM12.annotation.gff3.gz'
        mapping_url = 'data-raw/gencode.vM12.metadata.EntrezGene.gz'
        organism = 'Mus musculus'
        txdb_version = get_package_version('TxDb.Mmusculus.UCSC.mm10.knownGene')
        orgdb_version = get_package_version('org.Mm.eg.db')
		gencode_version = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz'
		mapping_version = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.metadata.EntrezGene.gz'
    } else if (genome == 'rn4') {
        # Gives ENSEMBL IDs
        require(TxDb.Rnorvegicus.UCSC.rn4.ensGene)
        require(org.Rn.eg.db)
        txdb = TxDb.Rnorvegicus.UCSC.rn4.ensGene
        egSYMBOL = org.Rn.egSYMBOL
        egENSEMBL2EG = org.Rn.egENSEMBL2EG
        egGENENAME = org.Rn.egGENENAME
        gencode_url = NULL
        mapping_url = NULL
        organism = 'Rattus norvegicus'
        txdb_version = get_package_version('TxDb.Rnorvegicus.UCSC.rn4.ensGene')
        orgdb_version = get_package_version('org.Rn.eg.db')
		gencode_version = NULL
		mapping_version = NULL
    } else if (genome == 'rn5') {
        # Gives Entrez IDs
        require(TxDb.Rnorvegicus.UCSC.rn5.refGene)
        require(org.Rn.eg.db)
        txdb = TxDb.Rnorvegicus.UCSC.rn5.refGene
        egSYMBOL = org.Rn.egSYMBOL
        egENSEMBL2EG = NULL
        egGENENAME = org.Rn.egGENENAME
        gencode_url = NULL
        mapping_url = NULL
        organism = 'Rattus norvegicus'
        txdb_version = get_package_version('TxDb.Rnorvegicus.UCSC.rn5.refGene')
        orgdb_version = get_package_version('org.Rn.eg.db')
		gencode_version = NULL
		mapping_version = NULL
    } else if (genome == 'rn6') {
        # Gives Entrez IDs
        require(TxDb.Rnorvegicus.UCSC.rn6.refGene)
        require(org.Rn.eg.db)
        txdb = TxDb.Rnorvegicus.UCSC.rn6.refGene
        egSYMBOL = org.Rn.egSYMBOL
        egENSEMBL2EG = NULL
        egGENENAME = org.Rn.egGENENAME
        gencode_url = NULL
        mapping_url = NULL
        organism = 'Rattus norvegicus'
        txdb_version = get_package_version('TxDb.Rnorvegicus.UCSC.rn6.refGene')
        orgdb_version = get_package_version('org.Rn.eg.db')
		gencode_version = NULL
		mapping_version = NULL
    } else if (genome == 'dm3') {
        # Gives ENSEMBL IDs
        require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
        require(org.Dm.eg.db)
        txdb = TxDb.Dmelanogaster.UCSC.dm3.ensGene
        egSYMBOL = org.Dm.egSYMBOL
        egENSEMBL2EG = org.Dm.egENSEMBL2EG
        egGENENAME = org.Dm.egGENENAME
        gencode_url = NULL
        mapping_url = NULL
        organism = 'Drosophila melanogaster'
        txdb_version = get_package_version('TxDb.Dmelanogaster.UCSC.dm3.ensGene')
        orgdb_version = get_package_version('org.Dm.eg.db')
		gencode_version = NULL
		mapping_version = NULL
    } else if (genome == 'dm6') {
        # Gives ENSEMBL IDs
        require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
        require(org.Dm.eg.db)
        txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene
        egSYMBOL = org.Dm.egSYMBOL
        egENSEMBL2EG = org.Dm.egENSEMBL2EG
        egGENENAME = org.Dm.egGENENAME
        gencode_url = NULL
        mapping_url = NULL
        organism = 'Drosophila melanogaster'
        txdb_version = get_package_version('TxDb.Dmelanogaster.UCSC.dm6.ensGene')
        orgdb_version = get_package_version('org.Dm.eg.db')
		gencode_version = NULL
		mapping_version = NULL
    } else {
        stop('Select a valid genome')
    }

    message(sprintf('BUILDING locus definitions for %s', genome))

    if(!is.null(gencode_url)) {
        ### GENCODE annotations
        # 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz'
        gencode = rtracklayer::readGFF(gencode_url)

        ### ENSEMBL transcript ID to Entrez ID mapping
        # 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz'
        ens2eg = readr::read_tsv(mapping_url, col_names = c('ens_id','gene_id'))
    } else {
        gencode = NULL
        ens2eg = NULL
    }

    if(is.null(egENSEMBL2EG)) {
        ### Build Entrez ID to gene symbol mapping
            # Get the gene symbol that are mapped to an entrez gene identifiers
            mapped_genes = AnnotationDbi::mappedkeys(egSYMBOL)
            # Convert to a data.frame
            eg2symbol = as.data.frame(egSYMBOL[mapped_genes])
            eg2symbol$gene_id = as.integer(eg2symbol$gene_id)

        ### Make a null ensembl2eg
            ensembl2eg = NULL
    } else {
        ### Build ENSEMBL ID to Entrez ID mapping
            mapped_genes = AnnotationDbi::mappedkeys(egENSEMBL2EG)
            ensembl2eg = as.data.frame(egENSEMBL2EG[mapped_genes])
            ensembl2eg$gene_id = as.integer(ensembl2eg$gene_id)
        ### Build Entrez ID table to gene symbol mapping
            mapped_genes = AnnotationDbi::mappedkeys(egSYMBOL)
            eg2symbol = as.data.frame(egSYMBOL[mapped_genes])
            eg2symbol$gene_id = as.integer(eg2symbol$gene_id)
    }

    ### Filter txdb for canonical chromosomes
        seqs = seqlevels(txdb)
        seqs = seqs[!(grepl('gl',seqs) | grepl('hap',seqs) | grepl('alt',seqs) | grepl('random',seqs) | grepl('chrUn',seqs))]
        seqlevels(txdb) = seqs

    ################################################################################
    ### Build up gr_transcripts for nearest_tss, nearest_gene, and all nkb definitions
    ################################################################################

        gr_transcripts = build_transcript_gr(txdb, eg2symbol, ensembl2eg)
        gr_transcripts = filter_gencode_types(gr = gr_transcripts, gencode = gencode, ens2eg = ens2eg)
        gr_transcripts = filter_readthrough_transcripts(gr = gr_transcripts, egGENENAME = egGENENAME)

        gr_transcripts = add_ntss_precursors(gr_transcripts)
        gr_transcripts = add_ngene_precursors(gr_transcripts)

        gr_transcripts = add_nkb_precursors(gr_transcripts, width = 1000)
        gr_transcripts = add_nkb_precursors(gr_transcripts, width = 5000)
        gr_transcripts = add_nkb_precursors(gr_transcripts, width = 10000)

    ################################################################################
    ### Build up gr_exons and gr_introns for exon and intron definitions
    ################################################################################

        gr_exons = build_exons_introns_gr(txdb, eg2symbol, ensembl2eg, type='exons')
        gr_exons = filter_gencode_types(gr = gr_exons, gencode = gencode, ens2eg = ens2eg)
        gr_exons = filter_readthrough_transcripts(gr = gr_exons, egGENENAME = egGENENAME)

        gr_exons = reduce_gr(gr_exons, eg2symbol)


        gr_introns = build_exons_introns_gr(txdb, eg2symbol, ensembl2eg, type='introns')
        gr_introns = filter_gencode_types(gr = gr_introns, gencode = gencode, ens2eg = ens2eg)
        gr_introns = filter_readthrough_transcripts(gr = gr_introns, egGENENAME = egGENENAME)

        gr_introns = reduce_gr(gr_introns, eg2symbol)

    ################################################################################
    ### Examine coverage of exons and introns
    ################################################################################

        cov_exons = GenomicRanges::coverage(gr_exons)
        all_cov_exons = c()
        for(chr in names(cov_exons)) {
            all_cov_exons = c(all_cov_exons, cov_exons[[chr]]@values)
        }
        table(all_cov_exons)

        cov_introns = GenomicRanges::coverage(gr_introns)
        all_cov_introns = c()
        for(chr in names(cov_introns)) {
            all_cov_introns = c(all_cov_introns, cov_introns[[chr]]@values)
        }
        table(all_cov_introns)

    ################################################################################
    ### Build locus definitions
    ################################################################################

    ### nearest_tss
        message('Constructing nearest_tss ldef...')
        ldef_ntss = build_ldef_nearest_tss(gr_transcripts, eg2symbol)

    ### nearest_gene
        message('Constructing nearest_gene ldef...')
        ldef_ngene = build_ldef_nearest_gene(gr_transcripts, eg2symbol)

    ### exons
        message('Constructing exon ldef...')
        ldef_exons = gr_exons

    ### introns
        message('Constructing intron ldef...')
        ldef_introns = gr_introns

    ### 1kb related definitions
        message('Constructing <1kb ldef...')
        ldef_1kb = build_ldef_nkb(gr_transcripts, width = 1000, context = 'inside', eg2symbol = eg2symbol)
        message('Constructing >1kb upstream ldef...')
        ldef_1kb_upstream = build_ldef_nkb(gr_transcripts, width = 1000, context = 'upstream', eg2symbol = eg2symbol)
        message('Constructing >1kb ldef...')
        ldef_1kb_outside = build_ldef_nkb(gr_transcripts, width = 1000, context = 'outside', eg2symbol = eg2symbol)

    ### 5kb related definitions
        message('Constructing <5kb ldef...')
        ldef_5kb = build_ldef_nkb(gr_transcripts, width = 5000, context = 'inside', eg2symbol = eg2symbol)
            message('Constructing >5kb upstream ldef...')
        ldef_5kb_upstream = build_ldef_nkb(gr_transcripts, width = 5000, context = 'upstream', eg2symbol = eg2symbol)
            message('Constructing >5kb ldef...')
        ldef_5kb_outside = build_ldef_nkb(gr_transcripts, width = 5000, context = 'outside', eg2symbol = eg2symbol)

    ### 10kb related definitions
        message('Constructing <10kb ldef...')
        ldef_10kb = build_ldef_nkb(gr_transcripts, width = 10000, context = 'inside', eg2symbol = eg2symbol)
            message('Constructing >10kb upstream ldef...')
        ldef_10kb_upstream = build_ldef_nkb(gr_transcripts, width = 10000, context = 'upstream', eg2symbol = eg2symbol)
            message('Constructing >10kb ldef...')
        ldef_10kb_outside = build_ldef_nkb(gr_transcripts, width = 10000, context = 'outside', eg2symbol = eg2symbol)

    ################################################################################
    ### Construct tss.[genome_build].RData
    ################################################################################

        to_tss_rdata(gr = gr_transcripts, genome = genome, txdb_version, orgdb_version, gencode_version, mapping_version)

    ################################################################################
    ### Construct LocusDefinition class objects for each
    ################################################################################

        ldef_gr_to_LocusDefinition(ldef_gr = ldef_ntss, genome = genome, organism = organism, ldef_name = 'nearest_tss', txdb_version, orgdb_version, gencode_version, mapping_version)
        ldef_gr_to_LocusDefinition(ldef_gr = ldef_ngene, genome = genome, organism = organism, ldef_name = 'nearest_gene', txdb_version, orgdb_version, gencode_version, mapping_version)

        ldef_gr_to_LocusDefinition(ldef_gr = ldef_exons, genome = genome, organism = organism, ldef_name = 'exon', txdb_version, orgdb_version, gencode_version, mapping_version)
        ldef_gr_to_LocusDefinition(ldef_gr = ldef_introns, genome = genome, organism = organism, ldef_name = 'intron', txdb_version, orgdb_version, gencode_version, mapping_version)

        ldef_gr_to_LocusDefinition(ldef_gr = ldef_1kb, genome = genome, organism = organism, ldef_name = '1kb', txdb_version, orgdb_version, gencode_version, mapping_version)
        ldef_gr_to_LocusDefinition(ldef_gr = ldef_1kb_upstream, genome = genome, organism = organism, ldef_name = '1kb_outside_upstream', txdb_version, orgdb_version, gencode_version, mapping_version)
        ldef_gr_to_LocusDefinition(ldef_gr = ldef_1kb_outside, genome = genome, organism = organism, ldef_name = '1kb_outside', txdb_version, orgdb_version, gencode_version, mapping_version)

        ldef_gr_to_LocusDefinition(ldef_gr = ldef_5kb, genome = genome, organism = organism, ldef_name = '5kb', txdb_version, orgdb_version, gencode_version, mapping_version)
        ldef_gr_to_LocusDefinition(ldef_gr = ldef_5kb_upstream, genome = genome, organism = organism, ldef_name = '5kb_outside_upstream', txdb_version, orgdb_version, gencode_version, mapping_version)
        ldef_gr_to_LocusDefinition(ldef_gr = ldef_5kb_outside, genome = genome, organism = organism, ldef_name = '5kb_outside', txdb_version, orgdb_version, gencode_version, mapping_version)

        ldef_gr_to_LocusDefinition(ldef_gr = ldef_10kb, genome = genome, organism = organism, ldef_name = '10kb', txdb_version, orgdb_version, gencode_version, mapping_version)
        ldef_gr_to_LocusDefinition(ldef_gr = ldef_10kb_upstream, genome = genome, organism = organism, ldef_name = '10kb_outside_upstream', txdb_version, orgdb_version, gencode_version, mapping_version)
        ldef_gr_to_LocusDefinition(ldef_gr = ldef_10kb_outside, genome = genome, organism = organism, ldef_name = '10kb_outside', txdb_version, orgdb_version, gencode_version, mapping_version)

    ################################################################################
    ### Output for sanity checks in Genome Browser
    ################################################################################

    ### ldefs to bed files for checking
        # ldef_gr_to_bed_file(ldef_ntss, sprintf('%s_nearest_tss',genome))
        # ldef_gr_to_bed_file(ldef_ngene, sprintf('%s_nearest_gene',genome))
        #
        # ldef_gr_to_bed_file(ldef_exons, sprintf('%s_exon',genome))
        # ldef_gr_to_bed_file(ldef_introns,sprintf('%s_intron',genome))
        #
        # ldef_gr_to_bed_file(ldef_1kb, sprintf('%s_1kb',genome))
        # ldef_gr_to_bed_file(ldef_1kb_upstream, sprintf('%s_1kb_upstream',genome))
        # ldef_gr_to_bed_file(ldef_1kb_outside, sprintf('%s_1kb_outside',genome))
        #
        # ldef_gr_to_bed_file(ldef_5kb, sprintf('%s_5kb',genome))
        # ldef_gr_to_bed_file(ldef_5kb_upstream, sprintf('%s_5kb_upstream',genome))
        # ldef_gr_to_bed_file(ldef_5kb_outside, sprintf('%s_5kb_outside',genome))
        #
        # ldef_gr_to_bed_file(ldef_10kb, sprintf('%s_10kb',genome))
        # ldef_gr_to_bed_file(ldef_10kb_upstream, sprintf('%s_10kb_upstream',genome))
        # ldef_gr_to_bed_file(ldef_10kb_outside, sprintf('%s_10kb_outside',genome))

    return(NULL)
}
