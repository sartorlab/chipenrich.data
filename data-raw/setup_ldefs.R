# The precede() and follow() functions are quiet weird so it is recommended that
# you look at a simple example and try to figure out what's going on.
#
# In particular, precede() and follow() return indices in subject for which the
# element of x is the direct precedent or follower, respectively, of. It would
# appear that precede() and follow() get it backwards.
#
# a = GRanges(
#     seqnames = c('chr1','chr1','chr1'),
#     ranges = IRanges(start = c(1,3,5), end = c(1,3,5)),
#     strand = '*')
#
# b = GRanges(
#     seqnames = c('chr1','chr1','chr1'),
#     ranges = IRanges(start = c(2,4,6), end = c(2,4,6)),
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

library(rtracklayer)
library(readr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

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
    build_transcript_gr = function(txdb, eg2symbol) {
        # Get transcripts with Entrez IDs and add gene symbols
        gr = transcripts(txdb, columns=c('gene_id'))
        gr$gene_id = as.integer(gr$gene_id)
        gr$symbol = eg2symbol[match(gr$gene_id, eg2symbol$gene_id), 'symbol']

        # Only admit transcripts with Entrez IDs
        gr = gr[!is.na(gr$gene_id)]

        gr = unique(gr)
        gr = sort(gr)

        return(gr)
    }

    ###
    ### Exons and introns
    ###
    build_exons_introns_gr = function(txdb, eg2symbol, type = c('exons','introns')) {
        type = match.arg(type)

        if(type == 'exons') {
            grl = exonsBy(txdb, by= 'tx', use.names = TRUE)
        } else {
            grl = intronsByTranscript(txdb, use.names = TRUE)
        }

        txname_rle = Rle(names(grl), elementNROWS(grl))
        txname_vec = as.character(txname_rle)

        gr = unlist(grl, use.names = FALSE)
        mcols(gr)$tx_name = txname_vec

        # Add Entrez ID and symbol
        # UCSC TXID and TXNAME to GENEID mapping (for introns and exons)
        id_maps = AnnotationDbi::select(txdb, keys = names(grl), columns = c('TXID','GENEID'), keytype = 'TXNAME')
        mcols(gr)$gene_id = id_maps[match(mcols(gr)$tx_name, id_maps$TXNAME), 'GENEID']
        mcols(gr)$symbol = eg2symbol[match(mcols(gr)$gene_id, eg2symbol$gene_id), 'symbol']
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

    filter_gencode_types = function(gr, gencode, ens2eg) {
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

        return(gr)
    }

################################################################################
### Filter readthrough symbols, e.g. APITD1-CORT
################################################################################

    filter_readthrough_transcripts = function(gr, egGENENAME) {
        ### Based on 'through' in gene name
            # Get a data.frame of Entrez IDs mapped to gene names
            mapped_genes = mappedkeys(egGENENAME)
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
### Use the seqinfo() of the subject GRanges to determine length of chromosomes.
### NOTE: seqinfo(x) should be the same as seqinfo(subject).
################################################################################

    subject_preceding_x = function(x, subject) {
        # Create GRanges of starts of chromosomes
        gr_chr_start = GRanges(
            seqnames = as.character(seqnames(seqinfo(subject))),
            ranges = IRanges(start = 1, end = 1),
            strand = '*',
            seqinfo = seqinfo(subject))

        # Add start ranges to subject so there are no NAs in follow()
        subject = sort(c(subject, gr_chr_start))

        # Determine indices of ranges in subject that precede items in x
        precede_idx = follow(x, subject, ignore.strand = TRUE)

        # Create initial draft of GRanges to return
        gr = granges(subject[precede_idx])

        return(gr)
    }

    subject_following_x = function(x, subject) {
        # Create GRanges of starts of chromosomes
        lengths = seqlengths(seqinfo(subject))
        gr_chr_end = GRanges(
            seqnames = as.character(seqnames(seqinfo(subject))),
            ranges = IRanges(start = lengths, end = lengths),
            strand = '*',
            seqinfo = seqinfo(subject))

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
    add_ntss_precursors = function(gr) {
        pos_idx = which(strand(gr) == '+')
        neg_idx = which(strand(gr) == '-')

        # Add an integer column for the TSS
        mcols(gr)$tss = 0
        mcols(gr[pos_idx])$tss = start(gr)[pos_idx]
        mcols(gr[neg_idx])$tss = end(gr)[neg_idx]

        # Add a GRanges column for the TSS
        mcols(gr)$gr_tss = GRanges(
            seqnames = seqnames(gr),
            ranges = IRanges(start = gr$tss, end = gr$tss),
            strand = strand(gr),
            seqinfo = seqinfo(gr))

        # Add GRanges columns for the following and preceding TSSs
        gr$tss_follow_tss = subject_following_x(gr$gr_tss, gr$gr_tss)
        gr$tss_precede_tss = subject_preceding_x(gr$gr_tss, gr$gr_tss)

        ### Determine boundaries of nearest_tss definitions. Essentially, midpoints btw TSSs

        # Add IRanges columns for starts and ends of locus def for transcripts
        # which will be used to take the midpoint
        gr$ntss_start_range = IRanges(start = start(gr$tss_precede_tss), end = start(gr$gr_tss))
        gr$ntss_end_range = IRanges(start = start(gr$gr_tss), end = start(gr$tss_follow_tss))

        # Add integer columns for start and end of locus def for transcript
        gr$ntss_start = mid(gr$ntss_start_range) + 1
        gr$ntss_end = mid(gr$ntss_end_range)

        # Add GRanges column of the final nearest_tss locus def
        gr$nearest_tss = GRanges(
            seqnames = seqnames(gr),
            ranges = IRanges(start = gr$ntss_start, end = gr$ntss_end),
            strand = '*')

        return(gr)
    }

################################################################################
### Add gr_left and gr_right mcols() to the input GRanges
################################################################################
    add_ngene_precursors = function(gr) {
        # Add GRanges for the left and right end of the feature
        mcols(gr)$gr_left = GRanges(
            seqnames = seqnames(gr),
            ranges = IRanges(start = start(gr), end = start(gr)),
            strand = strand(gr),
            seqinfo = seqinfo(gr))

        mcols(gr)$gr_right = GRanges(
            seqnames = seqnames(gr),
            ranges = IRanges(start = end(gr), end = end(gr)),
            strand = strand(gr),
            seqinfo = seqinfo(gr))

        # Combine the left and rights so it will find whichever is closest
        gr_combined = sort(c(gr$gr_right, gr$gr_left))

        # Add GRanges for precede and follow for left and right ends
        gr$precede_right = subject_preceding_x(gr$gr_right, gr_combined)
        gr$follow_right = subject_following_x(gr$gr_right, gr_combined)

        gr$precede_left = subject_preceding_x(gr$gr_left, gr_combined)
        gr$follow_left = subject_following_x(gr$gr_left, gr_combined)

        # Add IRanges with which to take midpoints for the nearest_gene definition
        ngene_left_start_range = IRanges(start = start(gr$precede_left), end = start(gr))
        ngene_left_end_range = IRanges(start = start(gr), end = start(gr$follow_left))

        ngene_right_start_range = IRanges(start = start(gr$precede_right), end = end(gr))
        ngene_right_end_range = IRanges(start = end(gr), end = start(gr$follow_right))

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

    add_nkb_precursors = function(gr, width) {
        ### Create shortcodes for nkbs
        if(width == 1000) {
            short_code = '1kb'
        } else if (width == 5000) {
            short_code = '5kb'
        } else if (width == 10000) {
            short_code = '10kb'
        }

        ### nkb base
            gr_nkb = flank(
                x = gr$gr_tss,
                width = width,
                both = TRUE,
                use.names = FALSE,
                ignore.strand = TRUE)
            seqinfo(gr_nkb) = seqinfo(gr_transcripts)

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
            seqinfo(gr_nkb_upstream) = seqinfo(gr)

            gr_nkb_downstream = c(gr_nkb_downstream_pos, gr_nkb_downstream_neg)
            seqinfo(gr_nkb_downstream) = seqinfo(gr)

        ### Add gr_nkb_fix, gr_nkb_upstream, and gr_nkb_downstream to mcols(gr)
            gr$gr_nkb = gr_nkb_fix
            gr$gr_nkb_upstream = gr_nkb_upstream
            gr$gr_nkb_downstream = gr_nkb_downstream

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

    reduce_gr = function(gr, eg2symbol) {
        message('Splitting on gene_id into GRangesList...')
        tmp_grl = splitAsList(gr, gr$gene_id)

        message('Reducing within gene_id...')
        tmp_grl = reduce(tmp_grl)

        # Get the correct vector of the gene_id from the GRangesList
        gene_id_rle = S4Vectors::Rle(names(tmp_grl), S4Vectors::elementNROWS(tmp_grl))
        gene_id_vec = as.integer(gene_id_rle)

        message('Reconstructing GRanges...')
        gr = unlist(tmp_grl, use.names = FALSE)

        # Add back the gene_ids and symbols
        message('Appending gene_id and symbol...')
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

    build_ldef_nearest_tss = function(gr, eg2symbol) {
        # Extract ntss
        ldef = gr$nearest_tss

        # Add gene_id and symbols
        ldef$gene_id = gr$gene_id
        ldef$symbol = gr$symbol

        # Inherit seqinfo()
        seqinfo(ldef) = seqinfo(gr)

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

    build_ldef_nearest_gene = function(gr, eg2symbol) {
        ### Create the nearest_gene definition by concatenating the GRanges consisting of:
        ### [gr$precede_left, start(gr)] and [start(gr), gr$follow_left]
        ### [gr$precede_right, end(gr)] and [end(gr), gr$follow_right]
        ### In other words, we are expanding around each transcript's left and right ends,
        ### ranges up until the preceding and following. There will initially be some
        ### redundancy that will be fixed with reduce() in a GRangesList on the gene_id.

        # Build the nearest_gene definition
        ldef = c(
            GRanges(
                seqnames = seqnames(gr),
                ranges = IRanges(start = gr$ngene_left_start, end = gr$ngene_left_end),
                strand = '*',
                gene_id = gr$gene_id,
                symbol = gr$symbol),
            GRanges(
                seqnames = seqnames(gr),
                ranges = IRanges(start = gr$ngene_right_start, end = gr$ngene_right_end),
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

    # @param width Is the width of the nkb definition, e.g. 1000, 5000, or 10000
    # @param context Is one of inside, upstream, or outside
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

    ldef_gr_to_bed_file = function(ldef_gr, ldef_name) {
        df = data.frame(ldef_gr, stringsAsFactors=F)
        df$strand = '.'
        df$score = 1000
        df$name = paste(df$gene_id, df$symbol, sep=':')

        df = df[,c('seqnames','start','end','name','score','strand')]
        write.table(df, file = sprintf('data/%s.bed', ldef_name), sep='\t', quote = F, col.names=F, row.names=F)
    }

    old_ldef_to_bed_file = function(old_ldef, ldef_name) {
        df = data.frame(old_ldef@granges, stringsAsFactors=F)
        df$strand = '.'
        df$score = 1000

        df = df[,c('seqnames','start','end','names','score','strand')]
        write.table(df, file = sprintf('data/%s.bed', ldef_name), sep='\t', quote = F, col.names=F, row.names=F)
    }

    ldef_gr_to_LocusDefinition = function(ldef_gr, genome, organism, ldef_name) {
        object = new("LocusDefinition")

        object@granges = ldef_gr
        object@dframe = data.frame(ldef_gr)
        object@genome.build = genome
        object@organism = organism

        save(object, file = sprintf('data/locusdef.%s.%s.RData', genome, ldef_name), compress = 'xz')
    }

################################################################################
### MAIN
################################################################################

### Establish data sources for base of definitions
    ### For the transcripts, TSSs, TESs, and Entrez Gene IDs
    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene

    ### Entrez ID to gene symbol mapping
    orgdb = org.Hs.egSYMBOL

### Establish data sources for filtering of definitions
    ### GENCODE data for filtering types of transcripts
    gencode_url = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz'
    mapping_url = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz'

    ### GENCODE annotations
    # 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gff3.gz'
    gencode = readGFF(gencode_url)

    ### ENSEMBL transcript ID to Entrez ID mapping
    # 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz'
    ens2eg = read_tsv(mapping_url, col_names = c('ens_id','gene_id'))

### Build Entrez ID to gene symbol mapping
    # Get the gene symbol that are mapped to an entrez gene identifiers
    mapped_genes = mappedkeys(orgdb)
    # Convert to a data.frame
    eg2symbol = as.data.frame(orgdb[mapped_genes])
    eg2symbol$gene_id = as.integer(eg2symbol$gene_id)

### Filter txdb for canonical chromosomes
    seqs = seqlevels(txdb)
    seqs = seqs[!(grepl('gl',seqs) | grepl('hap',seqs))]
    seqlevels(txdb) = seqs

################################################################################
### Build up gr_transcripts for nearest_tss, nearest_gene, and all nkb definitions
################################################################################

    gr_transcripts = build_transcript_gr(txdb, eg2symbol)
    gr_transcripts = filter_gencode_types(gr = gr_transcripts, gencode = gencode, ens2eg = ens2eg)
    gr_transcripts = filter_readthrough_transcripts(gr = gr_transcripts, egGENENAME = org.Hs.egGENENAME)

    gr_transcripts = add_ntss_precursors(gr_transcripts)
    gr_transcripts = add_ngene_precursors(gr_transcripts)

    gr_transcripts = add_nkb_precursors(gr_transcripts, width = 1000)
    gr_transcripts = add_nkb_precursors(gr_transcripts, width = 5000)
    gr_transcripts = add_nkb_precursors(gr_transcripts, width = 10000)

################################################################################
### Build up gr_exons and gr_introns for exon and intron definitions
################################################################################

    gr_exons = build_exons_introns_gr(txdb, eg2symbol, type='exons')
    gr_exons = filter_gencode_types(gr = gr_exons, gencode = gencode, ens2eg = ens2eg)
    gr_exons = filter_readthrough_transcripts(gr = gr_exons, egGENENAME = org.Hs.egGENENAME)

    gr_exons = reduce_gr(gr_exons, eg2symbol)


    gr_introns = build_exons_introns_gr(txdb, eg2symbol, type='introns')
    gr_introns = filter_gencode_types(gr = gr_introns, gencode = gencode, ens2eg = ens2eg)
    gr_introns = filter_readthrough_transcripts(gr = gr_introns, egGENENAME = org.Hs.egGENENAME)

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
### Construct LocusDefinition class objects for each
################################################################################

genome = 'hg19'
organism = 'Homo Sapiens'

    ldef_gr_to_LocusDefinition(ldef_gr = ldef_ntss, genome = genome, organism = organism, ldef_name = 'nearest_tss')
    ldef_gr_to_LocusDefinition(ldef_gr = ldef_ngene, genome = genome, organism = organism, ldef_name = 'nearest_gene')

    ldef_gr_to_LocusDefinition(ldef_gr = ldef_exons, genome = genome, organism = organism, ldef_name = 'exon')
    ldef_gr_to_LocusDefinition(ldef_gr = ldef_introns, genome = genome, organism = organism, ldef_name = 'intron')

    ldef_gr_to_LocusDefinition(ldef_gr = ldef_1kb, genome = genome, organism = organism, ldef_name = '1kb')
    ldef_gr_to_LocusDefinition(ldef_gr = ldef_1kb_upstream, genome = genome, organism = organism, ldef_name = '1kb_outside_upstream')
    ldef_gr_to_LocusDefinition(ldef_gr = ldef_1kb_outside, genome = genome, organism = organism, ldef_name = '1kb_outside')

    ldef_gr_to_LocusDefinition(ldef_gr = ldef_5kb, genome = genome, organism = organism, ldef_name = '5kb')
    ldef_gr_to_LocusDefinition(ldef_gr = ldef_5kb_upstream, genome = genome, organism = organism, ldef_name = '5kb_outside_upstream')
    ldef_gr_to_LocusDefinition(ldef_gr = ldef_5kb_outside, genome = genome, organism = organism, ldef_name = '5kb_outside')

    ldef_gr_to_LocusDefinition(ldef_gr = ldef_10kb, genome = genome, organism = organism, ldef_name = '10kb')
    ldef_gr_to_LocusDefinition(ldef_gr = ldef_10kb_upstream, genome = genome, organism = organism, ldef_name = '10kb_outside_upstream')
    ldef_gr_to_LocusDefinition(ldef_gr = ldef_10kb_outside, genome = genome, organism = organism, ldef_name = '10kb_outside')

################################################################################
### Output for sanity checks in Genome Browser
################################################################################

### ldefs to bed files for checking
    ldef_gr_to_bed_file(ldef_ntss, 'new_nearest_tss')
    ldef_gr_to_bed_file(ldef_ngene, 'new_nearest_gene')

    ldef_gr_to_bed_file(ldef_exons, 'new_exons')
    ldef_gr_to_bed_file(ldef_introns, 'new_introns')

    ldef_gr_to_bed_file(ldef_1kb, 'new_1kb')
    ldef_gr_to_bed_file(ldef_1kb_upstream, 'new_1kb_upstream')
    ldef_gr_to_bed_file(ldef_1kb_outside, 'new_1kb_outside')

    ldef_gr_to_bed_file(ldef_5kb, 'new_5kb')
    ldef_gr_to_bed_file(ldef_5kb_upstream, 'new_5kb_upstream')
    ldef_gr_to_bed_file(ldef_5kb_outside, 'new_5kb_outside')

    ldef_gr_to_bed_file(ldef_10kb, 'new_10kb')
    ldef_gr_to_bed_file(ldef_10kb_upstream, 'new_10kb_upstream')
    ldef_gr_to_bed_file(ldef_10kb_outside, 'new_10kb_outside')

### old ldefs to bed files for checking
    library(chipenrich.data)
    data(locusdef.hg19.nearest_tss, package = 'chipenrich.data')
    data(locusdef.hg19.nearest_gene, package = 'chipenrich.data')
    data(locusdef.hg19.exon, package = 'chipenrich.data')
    data(locusdef.hg19.intron, package = 'chipenrich.data')
    data(locusdef.hg19.1kb, package = 'chipenrich.data')
    data(locusdef.hg19.5kb, package = 'chipenrich.data')
    data(locusdef.hg19.10kb, package = 'chipenrich.data')
    data(locusdef.hg19.10kb_and_more_upstream, package = 'chipenrich.data')

    old_ldef_to_bed_file(locusdef.hg19.nearest_tss, 'ce_nearest_tss')
    old_ldef_to_bed_file(locusdef.hg19.nearest_gene, 'ce_nearest_gene')
    old_ldef_to_bed_file(locusdef.hg19.exon, 'ce_exons')
    old_ldef_to_bed_file(locusdef.hg19.intron, 'ce_introns')
    old_ldef_to_bed_file(locusdef.hg19.1kb, 'ce_1kb')
    old_ldef_to_bed_file(locusdef.hg19.5kb, 'ce_5kb')
    old_ldef_to_bed_file(locusdef.hg19.10kb, 'ce_10kb')
    old_ldef_to_bed_file(locusdef.hg19.10kb_and_more_upstream, 'ce_10kb_upstream')

################################################################################
