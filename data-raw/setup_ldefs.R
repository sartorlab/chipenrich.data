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

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

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

        # Determine indices of ranges in subject that precede items in x
        precede_idx = follow(x, subject, ignore.strand = TRUE)

        # Determine which are NA. These are the x ranges at the beginning of a chromosome
        na_idx = which(is.na(precede_idx))

        # Reset the NA indices to themselves
        precede_idx[na_idx] = match(subject[na_idx], subject)

        # Create initial draft of GRanges to return
        tmp_gr = granges(subject[precede_idx])

        # Determine seqnames that are NA. NOTE: There should be at most one per chromosome
        # since there can only be one first or last range per chromosome.
        subject_na_seqnames = as.character(seqnames(subject[na_idx]))
        gr_chr_start_seqnames = as.character(seqnames(gr_chr_start))

        # Create a GRanges to actually fix the na_idx
        tmp_fix_gr = gr_chr_start[match(subject_na_seqnames, gr_chr_start_seqnames)]

        # Apply the fix
        tmp_gr[na_idx] = tmp_fix_gr

        return(tmp_gr)
    }

    subject_following_x = function(x, subject) {
        # Create GRanges of starts of chromosomes
        lengths = seqlengths(seqinfo(subject))
        gr_chr_end = GRanges(
            seqnames = as.character(seqnames(seqinfo(subject))),
            ranges = IRanges(start = lengths, end = lengths),
            strand = '*',
            seqinfo = seqinfo(subject))

        # Determine indices of ranges in subject that precede items in x
        follow_idx = precede(x, subject, ignore.strand = TRUE)

        # Determine which are NA. These are the x ranges at the beginning of a chromosome
        na_idx = which(is.na(follow_idx))

        # Reset the NA indices to themselves
        follow_idx[na_idx] = match(subject[na_idx], subject)

        # Create initial draft of GRanges to return
        tmp_gr = granges(subject[follow_idx])

        # Determine seqnames that are NA. NOTE: There should be at most one per chromosome
        # since there can only be one first or last range per chromosome.
        subject_na_seqnames = as.character(seqnames(subject[na_idx]))
        gr_chr_end_seqnames = as.character(seqnames(gr_chr_end))

        # Create a GRanges to actually fix the na_idx
        tmp_fix_gr = gr_chr_end[match(subject_na_seqnames, gr_chr_end_seqnames)]

        # Apply the fix
        tmp_gr[na_idx] = tmp_fix_gr

        return(tmp_gr)
    }

    reduce_gr = function(gr) {
        message('Splitting on gene_id...')
        tmp_grl = splitAsList(gr, gr$gene_id)
        message('Reducing within gene_id...')
        tmp_grl = endoapply(tmp_grl, function(g){
            tmp = reduce(g)
            tmp$gene_id = unique(g$gene_id)
            tmp$symbol = unique(g$symbol)
            return(tmp)
        })
        message('Reconstructing GRanges...')
        tmp_gr = unlist(tmp_grl, use.names = FALSE)
        tmp_gr = sort(tmp_gr)
        return(tmp_gr)
    }

################################################################################
### Setup objects
################################################################################

###
### Establish databases
###

    # For the transcripts, TSSs, TESs, and Entrez Gene IDs
    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene

    ### Entrez ID to gene symbol mapping
    orgdb = org.Hs.egSYMBOL

###
### Filter for canonical chromosomes
###

    seqs = seqlevels(txdb)
    seqs = seqs[!(grepl('gl',seqs) | grepl('hap',seqs))]
    seqlevels(txdb) = seqs

###
### Build GRanges of transcripts, remove those with NA gene_id and make unique
### according to location + gene_id.
###

    # Build the transcripts
    gr = transcripts(txdb, columns=c('gene_id'))

    # Force the gene_id to integer from CharacterList
    # NOTE: Worth a check that mean(sapply(gr$gene_id, length)) = 1
    # so we're sure a transcript doesn't have multiple gene_ids.
    gr$gene_id = as.integer(gr$gene_id)

    # Force the range to have a gene_id
    # Can't do chipenrich without it
    gr = gr[!is.na(gr$gene_id)]

    # Enforce uniqueness on location + gene_id
    gr = unique(gr)

###
### Attach gene symbols
###

    # Get the gene symbol that are mapped to an entrez gene identifiers
    mapped_genes = mappedkeys(orgdb)
    # Convert to a data.frame
    eg2symbol = as.data.frame(orgdb[mapped_genes])
    eg2symbol$gene_id = as.integer(eg2symbol$gene_id)

    # Join eg2symbol on the gene_id column (allow NAs)
    gr$symbol = eg2symbol[match(gr$gene_id, eg2symbol$gene_id), 'symbol']

###
### Establish correct TSS and TES for each transcript according to strand.
### i.e. TSS = (start && +) && (end && -) and TES = (end && -) && (start && +)
### NOTE: Important for all the definitions.
###

    pos_idx = which(strand(gr) == '+')
    neg_idx = which(strand(gr) == '-')

    mcols(gr)$tss = 0
    mcols(gr[pos_idx])$tss = start(gr)[pos_idx]
    mcols(gr[neg_idx])$tss = end(gr)[neg_idx]

    mcols(gr)$tes = 0
    mcols(gr[pos_idx])$tes = end(gr)[pos_idx]
    mcols(gr[neg_idx])$tes = start(gr)[neg_idx]

    mcols(gr)$gr_tss = GRanges(
        seqnames = seqnames(gr),
        ranges = IRanges(start = gr$tss, end = gr$tss),
        strand = strand(gr),
        seqinfo = seqinfo(gr))

    mcols(gr)$gr_tes = GRanges(
        seqnames = seqnames(gr),
        ranges = IRanges(start = gr$tes, end = gr$tes),
        strand = strand(gr),
        seqinfo = seqinfo(gr))

###
### Establish left and right boundaries for each transcript
### NOTE: Important for nearets_gene definition.
###

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

###
### Establish nkb flanks around TSSs for each transcript
###

    # 1kb
    mcols(gr)$onekb = flank(
        x = gr$gr_tss,
        width = 1000,
        both = TRUE,
        use.names = FALSE,
        ignore.strand = TRUE)

    # 5kb
    mcols(gr)$fivekb = flank(
        x = gr$gr_tss,
        width = 5000,
        both = TRUE,
        use.names = FALSE,
        ignore.strand = TRUE)

    # 10kb
    mcols(gr)$tenkb = flank(
        x = gr$gr_tss,
        width = 10000,
        both = TRUE,
        use.names = FALSE,
        ignore.strand = TRUE)

###
### Filter for snoRNA, psueudogenes, lncRNA (?)
###

    # Look for snoRNA/miRNAs (COME BACK TO THIS)

    # snoquery = query(ah, c('sno/miRNA', 'Homo sapiens'))
    # snoid = names(snoquery)[snoquery$genome == 'hg19']
    # snodb = ah[[snoid]]
    #
    # # Match snoRNAs/miRNAs by location
    # found_snomi = which(match(gr, snodb) != 'NA')
    # found_snomi_idx = found_snomi[!is.na(found_snomi)]

################################################################################
### Build intermediate boundaries as GRanges to construct locus definitions
################################################################################

### In all that follows X is the object of interest and v is the object chosen
### by the procedure called. [ and ] denote TSSs and TESs and - is intergenic space.

################################################################################
### Determine following and preceding TSSs relative to all TSSs

###
### Determine TSSs following TSSs
###                        X                 v
### ------tss-------------TSS---------------tss------tss--------tsss
    gr$tss_follow_tss = subject_following_x(gr$gr_tss, gr$gr_tss)

###
### Determine TSSs preceding TSSs
###        v               X
### ------tss-------------TSS---------------tss------tss--------tsss
    gr$tss_precede_tss = subject_preceding_x(gr$gr_tss, gr$gr_tss)

################################################################################
### Determine following and preceding left and right ends of transcripts for
### right ends of transcripts.

###
### Determine right preceding right
###            v                       X
### -----[     ]--------------[        ]------------[       ]------------
###
###                                v    X
### -----[     ]--------[     [    ]    ]------------[       ]------------
    gr$right_precede_right = subject_preceding_x(gr$gr_right, gr$gr_right)

    # Get the distance
    gr$dist_right_precede_right = distance(gr$gr_right, gr$right_precede_right, ignore.strand = TRUE)

###
### Determine left preceding right
###            v                       X
### -----[     ]--------------[        ]------------[       ]------------
###
###                           v         X
### -----[     ]--------[     [    ]    ]------------[       ]------------
    gr$left_precede_right = subject_preceding_x(gr$gr_right, gr$gr_left)

    # Get the distance
    gr$dist_left_precede_right = distance(gr$gr_right, gr$left_precede_right, ignore.strand = TRUE)

###
### Decide whether preceding left or preceding right precedes right based on minimum distance
###                           v
###            pr             pl       X
### -----[     ]--------------[        ]------------[       ]------------
###                                v
###                           pl   pr   X
### -----[     ]--------[     [    ]    ]------------[       ]------------

    # Start by making the right the winner
    gr$precede_right = gr$right_precede_right

    # Get the indices where left is closer than right
    precede_right_left_closer_idx = which(gr$dist_right_precede_right > gr$dist_left_precede_right)
    # Replace gr$precede_right with the left where the preceding left is closer than the preceding right
    gr$precede_right[precede_right_left_closer_idx] = gr$left_precede_right[precede_right_left_closer_idx]

######
### NOTE: Cartoons for the other combinations are ommitted, but you get the idea.

###
### Determine right following right
###
    gr$right_follow_right = subject_following_x(gr$gr_right, gr$gr_right)

    # Get the distance
    gr$dist_right_follow_right = distance(gr$gr_right, gr$right_follow_right, ignore.strand = TRUE)

###
### Determine left following right
###
    gr$left_follow_right = subject_following_x(gr$gr_right, gr$gr_left)

    # Get the distance
    gr$dist_left_follow_right = distance(gr$gr_right, gr$left_follow_right, ignore.strand = TRUE)

###
### Decide whether following right or following left follows right based on minimum distance
###
    # Start by making the right the winner
    gr$follow_right = gr$right_follow_right

    # Get the indices where left is closer than right
    follow_right_left_closer_idx = which(gr$dist_right_follow_right > gr$dist_left_follow_right)
    # Replace gr$follow_right with the left where the following left is closer than the following right
    gr$follow_right[follow_right_left_closer_idx] = gr$left_follow_right[follow_right_left_closer_idx]

################################################################################
### Determine following and preceding left and right ends of transcripts for
### left ends of transcripts.

###
### Determine right preceding left
###
    gr$right_precede_left = subject_preceding_x(gr$gr_left, gr$gr_right)

    # Get the distance
    gr$dist_right_precede_left = distance(gr$gr_left, gr$right_precede_left, ignore.strand = TRUE)

###
### Determine left preceding left
###
    gr$left_precede_left = subject_preceding_x(gr$gr_left, gr$gr_left)

    # Get the distance
    gr$dist_left_precede_left = distance(gr$gr_left, gr$left_precede_left, ignore.strand = TRUE)

###
### Decide whether preceding right or preceding left precedes left based on minimum distance
###
    # Start by making the right the winner
    gr$precede_left = gr$right_precede_left

    # Get the indices where left is closer than right
    precede_left_left_closer_idx = which(gr$dist_right_precede_left > gr$dist_left_precede_left)
    # Replace gr$precede_left with the left where the preceding left is closer than the preceding right
    gr$precede_left[precede_left_left_closer_idx] = gr$left_precede_left[precede_left_left_closer_idx]

######

###
### Determine right following left
###
    gr$right_follow_left = subject_following_x(gr$gr_left, gr$gr_right)

    # Get the distance
    gr$dist_right_follow_left = distance(gr$gr_left, gr$right_follow_left, ignore.strand = TRUE)

###
### Determine left following left
###
    gr$left_follow_left = subject_following_x(gr$gr_left, gr$gr_left)

    # Get the distance
    gr$dist_left_follow_left = distance(gr$gr_left, gr$left_follow_left, ignore.strand = TRUE)

###
### Decide whether following right or following left follows left based on minimum distance
###
    # Start by making the right the winner
    gr$follow_left = gr$right_follow_left

    # Get the indices where left is closer than right
    follow_left_left_closer_idx = which(gr$dist_right_follow_left > gr$dist_left_follow_left)
    # Replace gr$follow_left with the left where the following left is closer than the following right
    gr$follow_left[follow_left_left_closer_idx] = gr$left_follow_left[follow_left_left_closer_idx]

################################################################################
### Determine boundaries of nearest_tss definitions. Essentially, midpoints btw TSSs

    # Setup IRanges from which to take the midpoints for the nearest_tss boundaries
    gr$ntss_start_range = IRanges(start = start(gr$tss_precede_tss), end = start(gr$gr_tss))
    gr$ntss_end_range = IRanges(start = start(gr$gr_tss), end = start(gr$tss_follow_tss))

    # The nearest_tss boundaries are the midpoints between nearest TSSs on either side
    # of the TSS of interest.
    gr$ntss_start = mid(gr$ntss_start_range) + 1
    gr$ntss_end = mid(gr$ntss_end_range)

################################################################################
### Build locus definitions
################################################################################

################################################################################
### Create the nearest_tss definition by using the midpoints of:
### [gr$tss_precede_tss, start(gr$gr_tss)] and [start(gr$gr_tss), gr$tss_follow_tss]
### as the start and endpoints of the nearest_tss range.

    # Build the nearest_tss definition
    gr_nearest_tss = GRanges(
        seqnames = seqnames(gr),
        ranges = IRanges(start = gr$ntss_start, end = gr$ntss_end),
        strand = '*',
        gene_id = gr$gene_id,
        symbol = gr$symbol)

    gr$nearest_tss = granges(gr_nearest_tss)

    # Enforce uniqueness of ranges/gene_id
    gr_nearest_tss = sort(gr_nearest_tss)
    gr_nearest_tss = unique(gr_nearest_tss)

################################################################################
### Create the nearest_gene definition by concatenating the GRanges consisting of:
### [gr$precede_left, start(gr)] and [start(gr), gr$follow_left]
### [gr$precede_right, end(gr)] and [end(gr), gr$follow_right]
### In other words, we are expanding around each transcript's left and right ends,
### ranges up until the preceding and following. There will initially be some
### redundancy that will be fixed with reduce() in a GRangesList on the gene_id.

    # Setup IRanges from which to take the midpoints for the nearest_gene boundaries
    gr$ngene_left_start_range = IRanges(start = start(gr$precede_left), end = start(gr))
    gr$ngene_left_end_range = IRanges(start = start(gr), end = start(gr$follow_left))

    gr$ngene_right_start_range = IRanges(start = start(gr$precede_right), end = end(gr))
    gr$ngene_right_end_range = IRanges(start = end(gr), end = start(gr$follow_right))

    # The nearest_gene boundaries are defined as the midpoints between the left
    # end of the range and the closest preceding TSS or TES and the right end of the
    # range and the closest following TSS or TES.
    gr$ngene_left_start = mid(gr$ngene_left_start_range) + 1
    gr$ngene_left_end = mid(gr$ngene_left_end_range)

    gr$ngene_right_start = mid(gr$ngene_right_start_range) + 1
    gr$ngene_right_end = mid(gr$ngene_right_end_range)

    # Build the nearest_gene definition
    gr_nearest_gene = c(
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
    gr_nearest_gene = sort(gr_nearest_gene)
    gr_nearest_gene = unique(gr_nearest_gene)

################################################################################

### Create 1kb definition by determining if the 1kb start or end exceeds
### ntss_start and ntss_end, respectively.
    onekb_start_less_ntss_start_idx = which(start(gr$onekb) < gr$ntss_start)
    onekb_end_greater_ntss_end_idx = which(end(gr$onekb) > gr$ntss_end)

    # Fix the 1kb ranges that go too far
    gr$onekb_fix = gr$onekb
    start(gr$onekb_fix)[onekb_start_less_ntss_start_idx] = gr$ntss_start[onekb_start_less_ntss_start_idx]
    end(gr$onekb_fix)[onekb_end_greater_ntss_end_idx] = gr$ntss_end[onekb_end_greater_ntss_end_idx]

    # Create the 1kb definition
    gr_onekb = gr$onekb_fix
    strand(gr_onekb) = '*'
    mcols(gr_onekb) = mcols(gr)[,c('gene_id','symbol')]

    gr_onekb = sort(gr_onekb)
    gr_onekb = unique(gr_onekb)

################################################################################

### Create 5kb definition by determining if the 5kb start or end exceeds
### ntss_start and ntss_end, respectively.
    fivekb_start_less_ntss_start_idx = which(start(gr$fivekb) < gr$ntss_start)
    fivekb_end_greater_ntss_end_idx = which(end(gr$fivekb) > gr$ntss_end)

    # Fix the 5kb ranges that go too far
    gr$fivekb_fix = gr$fivekb
    start(gr$fivekb_fix)[fivekb_start_less_ntss_start_idx] = gr$ntss_start[fivekb_start_less_ntss_start_idx]
    end(gr$fivekb_fix)[fivekb_end_greater_ntss_end_idx] = gr$ntss_end[fivekb_end_greater_ntss_end_idx]

    ### Create the 5kb definition
    gr_fivekb = gr$fivekb_fix
    strand(gr_fivekb) = '*'
    mcols(gr_fivekb) = mcols(gr)[,c('gene_id','symbol')]

    gr_fivekb = sort(gr_fivekb)
    gr_fivekb = unique(gr_fivekb)

################################################################################

### Create 10kb definition by determining if the 10kb start or end exceeds
### ntss_start and ntss_end, respectively.
    tenkb_start_less_ntss_start_idx = which(start(gr$tenkb) < gr$ntss_start)
    tenkb_end_greater_ntss_end_idx = which(end(gr$tenkb) > gr$ntss_end)

    # Fix the 10kb ranges that go too far
    gr$tenkb_fix = gr$tenkb
    start(gr$tenkb_fix)[tenkb_start_less_ntss_start_idx] = gr$ntss_start[tenkb_start_less_ntss_start_idx]
    end(gr$tenkb_fix)[tenkb_end_greater_ntss_end_idx] = gr$ntss_end[tenkb_end_greater_ntss_end_idx]

    ### Create a simple 10kb definition
    gr_tenkb = gr$tenkb_fix
    strand(gr_tenkb) = '*'
    mcols(gr_tenkb) = mcols(gr)[,c('gene_id','symbol')]

    gr_tenkb = sort(gr_tenkb)
    gr_tenkb = unique(gr_tenkb)

################################################################################
### Create the outside nkb definitions

###
### 1kb
###
    # Use restrict to get left and right side of 1kb def up to nearest_tss
    gr$outside_onekb_left = restrict(gr$nearest_tss, start = start(gr$nearest_tss), end = start(gr$onekb_fix))
    gr$outside_onekb_right = restrict(gr$nearest_tss, start = end(gr$onekb_fix), end = end(gr$nearest_tss))

    # Determine upstream ranges based on the strand
    gr_onekb_upstream_pos = gr$outside_onekb_left[strand(gr) == '+']
    mcols(gr_onekb_upstream_pos) = mcols(gr)[strand(gr) == '+', c('gene_id','symbol')]

    gr_onekb_upstream_neg = gr$outside_onekb_right[strand(gr) == '-']
    mcols(gr_onekb_upstream_neg) = mcols(gr)[strand(gr) == '-', c('gene_id','symbol')]

    # Determine downstream ranges based on the strand
    gr_onekb_downstream_pos = gr$outside_onekb_right[strand(gr) == '+']
    mcols(gr_onekb_downstream_pos) = mcols(gr)[strand(gr) == '+', c('gene_id','symbol')]

    gr_onekb_downstream_neg = gr$outside_onekb_left[strand(gr) == '-']
    mcols(gr_onekb_downstream_neg) = mcols(gr)[strand(gr) == '-', c('gene_id','symbol')]

    ### Create the >1kb upstream definition
    gr_onekb_upstream = c(gr_onekb_upstream_pos, gr_onekb_upstream_neg)
    gr_onekb_upstream = gr_onekb_upstream[which(width(gr_onekb_upstream) > 1)]

    gr_onekb_upstream = sort(gr_onekb_upstream)
    gr_onekb_upstream = unique(gr_onekb_upstream)

    ### Create the outside 1kb definition
    # Just need to create downstream and add it to the upstream
    gr_onekb_downstream = c(gr_onekb_downstream_pos, gr_onekb_downstream_neg)
    gr_onekb_downstream = gr_onekb_downstream[which(width(gr_onekb_downstream) > 1)]

    gr_onekb_outside = c(gr_onekb_upstream, gr_onekb_downstream)
    gr_onekb_outside = sort(gr_onekb_outside)
    gr_onekb_outside = unique(gr_onekb_outside)

###
### 5kb
###
    # Use restrict to get left and right side of 1kb def up to nearest_tss
    gr$outside_fivekb_left = restrict(gr$nearest_tss, start = start(gr$nearest_tss), end = start(gr$fivekb_fix))
    gr$outside_fivekb_right = restrict(gr$nearest_tss, start = end(gr$fivekb_fix), end = end(gr$nearest_tss))

    # Determine upstream ranges based on the strand
    gr_fivekb_upstream_pos = gr$outside_fivekb_left[strand(gr) == '+']
    mcols(gr_fivekb_upstream_pos) = mcols(gr)[strand(gr) == '+', c('gene_id','symbol')]

    gr_fivekb_upstream_neg = gr$outside_fivekb_right[strand(gr) == '-']
    mcols(gr_fivekb_upstream_neg) = mcols(gr)[strand(gr) == '-', c('gene_id','symbol')]

    # Determine downstream ranges based on the strand
    gr_fivekb_downstream_pos = gr$outside_fivekb_right[strand(gr) == '+']
    mcols(gr_fivekb_downstream_pos) = mcols(gr)[strand(gr) == '+', c('gene_id','symbol')]

    gr_fivekb_downstream_neg = gr$outside_fivekb_left[strand(gr) == '-']
    mcols(gr_fivekb_downstream_neg) = mcols(gr)[strand(gr) == '-', c('gene_id','symbol')]

    ### Create the >1kb upstream definition
    gr_fivekb_upstream = c(gr_fivekb_upstream_pos, gr_fivekb_upstream_neg)
    gr_fivekb_upstream = gr_fivekb_upstream[which(width(gr_fivekb_upstream) > 1)]

    gr_fivekb_upstream = sort(gr_fivekb_upstream)
    gr_fivekb_upstream = unique(gr_fivekb_upstream)

    ### Create the outside 1kb definition
    # Just need to create downstream and add it to the upstream
    gr_fivekb_downstream = c(gr_fivekb_downstream_pos, gr_fivekb_downstream_neg)
    gr_fivekb_downstream = gr_fivekb_downstream[which(width(gr_fivekb_downstream) > 1)]

    gr_fivekb_outside = c(gr_fivekb_upstream, gr_fivekb_downstream)
    gr_fivekb_outside = sort(gr_fivekb_outside)
    gr_fivekb_outside = unique(gr_fivekb_outside)

###
### 10kb
###
    # Use restrict to get left and right side of 1kb def up to nearest_tss
    gr$outside_tenkb_left = restrict(gr$nearest_tss, start = start(gr$nearest_tss), end = start(gr$tenkb_fix))
    gr$outside_tenkb_right = restrict(gr$nearest_tss, start = end(gr$tenkb_fix), end = end(gr$nearest_tss))

    # Determine upstream ranges based on the strand
    gr_tenkb_upstream_pos = gr$outside_tenkb_left[strand(gr) == '+']
    mcols(gr_tenkb_upstream_pos) = mcols(gr)[strand(gr) == '+', c('gene_id','symbol')]

    gr_tenkb_upstream_neg = gr$outside_tenkb_right[strand(gr) == '-']
    mcols(gr_tenkb_upstream_neg) = mcols(gr)[strand(gr) == '-', c('gene_id','symbol')]

    # Determine downstream ranges based on the strand
    gr_tenkb_downstream_pos = gr$outside_tenkb_right[strand(gr) == '+']
    mcols(gr_tenkb_downstream_pos) = mcols(gr)[strand(gr) == '+', c('gene_id','symbol')]

    gr_tenkb_downstream_neg = gr$outside_tenkb_left[strand(gr) == '-']
    mcols(gr_tenkb_downstream_neg) = mcols(gr)[strand(gr) == '-', c('gene_id','symbol')]

    ### Create the >1kb upstream definition
    gr_tenkb_upstream = c(gr_tenkb_upstream_pos, gr_tenkb_upstream_neg)
    gr_tenkb_upstream = gr_tenkb_upstream[which(width(gr_tenkb_upstream) > 1)]

    gr_tenkb_upstream = sort(gr_tenkb_upstream)
    gr_tenkb_upstream = unique(gr_tenkb_upstream)

    ### Create the outside 1kb definition
    # Just need to create downstream and add it to the upstream
    gr_tenkb_downstream = c(gr_tenkb_downstream_pos, gr_tenkb_downstream_neg)
    gr_tenkb_downstream = gr_tenkb_downstream[which(width(gr_tenkb_downstream) > 1)]

    gr_tenkb_outside = c(gr_tenkb_upstream, gr_tenkb_downstream)
    gr_tenkb_outside = sort(gr_tenkb_outside)
    gr_tenkb_outside = unique(gr_tenkb_outside)

################################################################################

### Reduce all definitions

gr_nearest_tss = reduce_gr(gr_nearest_tss)
gr_nearest_gene = reduce_gr(gr_nearest_gene)

gr_onekb = reduce_gr(gr_onekb)
gr_fivekb = reduce_gr(gr_fivekb)
gr_tenkb = reduce_gr(gr_tenkb)

gr_onekb_upstream = reduce_gr(gr_onekb_upstream)
gr_fivekb_upstream = reduce_gr(gr_fivekb_upstream)
gr_tenkb_upstream = reduce_gr(gr_tenkb_upstream)

gr_onekb_outside = reduce_gr(gr_onekb_outside)
gr_fivekb_outside = reduce_gr(gr_fivekb_outside)
gr_tenkb_outside = reduce_gr(gr_tenkb_outside)

################################################################################

### Output for sanity checks in Genome Browser

# ntss
df_ntss = data.frame(gr_nearest_tss, stringsAsFactors=F)
df_ntss$strand = '.'
df_ntss$score = 1000

df_ntss_symbol = df_ntss[,c('seqnames','start','end','symbol','score','strand')]
write.table(df_ntss_symbol, file = '~/Desktop/ntss_symbol.bed', sep='\t', quote = F, col.names=F, row.names=F)

df_ntss_geneid = df_ntss[,c('seqnames','start','end','gene_id','score','strand')]
write.table(df_ntss_geneid, file = '~/Desktop/ntss_geneid.bed', sep='\t', quote = F, col.names=F, row.names=F)

# ngene
df_ngene = data.frame(gr_nearest_gene, stringsAsFactors=F)
df_ngene$strand = '.'
df_ngene$score = 1000

df_ngene_symbol = df_ngene[,c('seqnames','start','end','symbol','score','strand')]
write.table(df_ngene_symbol, file = '~/Desktop/ngene_symbol.bed', sep='\t', quote = F, col.names=F, row.names=F)

df_ngene_geneid = df_ngene[,c('seqnames','start','end','gene_id','score','strand')]
write.table(df_ngene_geneid, file = '~/Desktop/ngene_geneid.bed', sep='\t', quote = F, col.names=F, row.names=F)

# onekb
df_onekb = data.frame(gr_onekb, stringsAsFactors=F)
df_onekb$strand = '.'
df_onekb$score = 1000

df_onekb_symbol = df_onekb[,c('seqnames','start','end','symbol','score','strand')]
write.table(df_onekb_symbol, file = '~/Desktop/onekb_symbol.bed', sep='\t', quote = F, col.names=F, row.names=F)

df_onekb_geneid = df_onekb[,c('seqnames','start','end','gene_id','score','strand')]
write.table(df_onekb_geneid, file = '~/Desktop/onekb_geneid.bed', sep='\t', quote = F, col.names=F, row.names=F)

# fivekb
df_fivekb = data.frame(gr_fivekb, stringsAsFactors=F)
df_fivekb$strand = '.'
df_fivekb$score = 1000

df_fivekb_symbol = df_fivekb[,c('seqnames','start','end','symbol','score','strand')]
write.table(df_fivekb_symbol, file = '~/Desktop/fivekb_symbol.bed', sep='\t', quote = F, col.names=F, row.names=F)

df_fivekb_geneid = df_fivekb[,c('seqnames','start','end','gene_id','score','strand')]
write.table(df_fivekb_geneid, file = '~/Desktop/fivekb_geneid.bed', sep='\t', quote = F, col.names=F, row.names=F)

# tenkb
df_tenkb = data.frame(gr_tenkb, stringsAsFactors=F)
df_tenkb$strand = '.'
df_tenkb$score = 1000

df_tenkb_symbol = df_tenkb[,c('seqnames','start','end','symbol','score','strand')]
write.table(df_tenkb_symbol, file = '~/Desktop/tenkb_symbol.bed', sep='\t', quote = F, col.names=F, row.names=F)

df_tenkb_geneid = df_tenkb[,c('seqnames','start','end','gene_id','score','strand')]
write.table(df_tenkb_geneid, file = '~/Desktop/tenkb_geneid.bed', sep='\t', quote = F, col.names=F, row.names=F)

# onekb_upstream
df_onekb_upstream = data.frame(gr_onekb_upstream, stringsAsFactors=F)
df_onekb_upstream$strand = '.'
df_onekb_upstream$score = 1000

df_onekb_upstream_symbol = df_onekb_upstream[,c('seqnames','start','end','symbol','score','strand')]
write.table(df_onekb_upstream_symbol, file = '~/Desktop/onekb_upstream_symbol.bed', sep='\t', quote = F, col.names=F, row.names=F)

df_onekb_upstream_geneid = df_onekb_upstream[,c('seqnames','start','end','gene_id','score','strand')]
write.table(df_onekb_upstream_geneid, file = '~/Desktop/onekb_upstream_geneid.bed', sep='\t', quote = F, col.names=F, row.names=F)

# fivekb_upstream
df_fivekb_upstream = data.frame(gr_fivekb_upstream, stringsAsFactors=F)
df_fivekb_upstream$strand = '.'
df_fivekb_upstream$score = 1000

df_fivekb_upstream_symbol = df_fivekb_upstream[,c('seqnames','start','end','symbol','score','strand')]
write.table(df_fivekb_upstream_symbol, file = '~/Desktop/fivekb_upstream_symbol.bed', sep='\t', quote = F, col.names=F, row.names=F)

df_fivekb_upstream_geneid = df_fivekb_upstream[,c('seqnames','start','end','gene_id','score','strand')]
write.table(df_fivekb_upstream_geneid, file = '~/Desktop/fivekb_upstream_geneid.bed', sep='\t', quote = F, col.names=F, row.names=F)

# tenkb_upstream
df_tenkb_upstream = data.frame(gr_tenkb_upstream, stringsAsFactors=F)
df_tenkb_upstream$strand = '.'
df_tenkb_upstream$score = 1000

df_tenkb_upstream_symbol = df_tenkb_upstream[,c('seqnames','start','end','symbol','score','strand')]
write.table(df_tenkb_upstream_symbol, file = '~/Desktop/tenkb_upstream_symbol.bed', sep='\t', quote = F, col.names=F, row.names=F)

df_tenkb_upstream_geneid = df_tenkb_upstream[,c('seqnames','start','end','gene_id','score','strand')]
write.table(df_tenkb_upstream_geneid, file = '~/Desktop/tenkb_upstream_geneid.bed', sep='\t', quote = F, col.names=F, row.names=F)

# onekb_outside
df_onekb_outside = data.frame(gr_onekb_outside, stringsAsFactors=F)
df_onekb_outside$strand = '.'
df_onekb_outside$score = 1000

df_onekb_outside_symbol = df_onekb_outside[,c('seqnames','start','end','symbol','score','strand')]
write.table(df_onekb_outside_symbol, file = '~/Desktop/onekb_outside_symbol.bed', sep='\t', quote = F, col.names=F, row.names=F)

df_onekb_outside_geneid = df_onekb_outside[,c('seqnames','start','end','gene_id','score','strand')]
write.table(df_onekb_outside_geneid, file = '~/Desktop/onekb_outside_geneid.bed', sep='\t', quote = F, col.names=F, row.names=F)

# fivekb_outside
df_fivekb_outside = data.frame(gr_fivekb_outside, stringsAsFactors=F)
df_fivekb_outside$strand = '.'
df_fivekb_outside$score = 1000

df_fivekb_outside_symbol = df_fivekb_outside[,c('seqnames','start','end','symbol','score','strand')]
write.table(df_fivekb_outside_symbol, file = '~/Desktop/fivekb_outside_symbol.bed', sep='\t', quote = F, col.names=F, row.names=F)

df_fivekb_outside_geneid = df_fivekb_outside[,c('seqnames','start','end','gene_id','score','strand')]
write.table(df_fivekb_outside_geneid, file = '~/Desktop/fivekb_outside_geneid.bed', sep='\t', quote = F, col.names=F, row.names=F)

# tenkb_outside
df_tenkb_outside = data.frame(gr_tenkb_outside, stringsAsFactors=F)
df_tenkb_outside$strand = '.'
df_tenkb_outside$score = 1000

df_tenkb_outside_symbol = df_tenkb_outside[,c('seqnames','start','end','symbol','score','strand')]
write.table(df_tenkb_outside_symbol, file = '~/Desktop/tenkb_outside_symbol.bed', sep='\t', quote = F, col.names=F, row.names=F)

df_tenkb_outside_geneid = df_tenkb_outside[,c('seqnames','start','end','gene_id','score','strand')]
write.table(df_tenkb_outside_geneid, file = '~/Desktop/tenkb_outside_geneid.bed', sep='\t', quote = F, col.names=F, row.names=F)

################################################################################

### Output old CE ldefs for sanity checks in Genome Browser

library(chipenrich.data)
data(locusdef.hg19.nearest_tss, package = 'chipenrich.data')
data(locusdef.hg19.nearest_gene, package = 'chipenrich.data')
data(locusdef.hg19.1kb, package = 'chipenrich.data')
data(locusdef.hg19.5kb, package = 'chipenrich.data')
data(locusdef.hg19.10kb, package = 'chipenrich.data')

ce_ntss = locusdef.hg19.nearest_tss@dframe
ce_ntss$strand = '.'
ce_ntss$score = 1000
ce_ntss = ce_ntss[,c('chrom','start','end','geneid','score','strand')]
write.table(ce_ntss, file = '~/Desktop/ce_ntss.bed', sep='\t', quote = F, col.names=F, row.names=F)

ce_ngene = locusdef.hg19.nearest_gene@dframe
ce_ngene$strand = '.'
ce_ngene$score = 1000
ce_ngene = ce_ngene[,c('chrom','start','end','geneid','score','strand')]
write.table(ce_ngene, file = '~/Desktop/ce_ngene.bed', sep='\t', quote = F, col.names=F, row.names=F)

ce_1kb = locusdef.hg19.1kb@dframe
ce_1kb$strand = '.'
ce_1kb$score = 1000
ce_1kb = ce_1kb[,c('chrom','start','end','geneid','score','strand')]
write.table(ce_1kb, file = '~/Desktop/ce_1kb.bed', sep='\t', quote = F, col.names=F, row.names=F)

ce_5kb = locusdef.hg19.5kb@dframe
ce_5kb$strand = '.'
ce_5kb$score = 1000
ce_5kb = ce_5kb[,c('chrom','start','end','geneid','score','strand')]
write.table(ce_5kb, file = '~/Desktop/ce_5kb.bed', sep='\t', quote = F, col.names=F, row.names=F)

ce_10kb = locusdef.hg19.10kb@dframe
ce_10kb$strand = '.'
ce_10kb$score = 1000
ce_10kb = ce_10kb[,c('chrom','start','end','geneid','score','strand')]
write.table(ce_10kb, file = '~/Desktop/ce_10kb.bed', sep='\t', quote = F, col.names=F, row.names=F)

ce_10kb_upstream = locusdef.hg19.10kb_and_more_upstream@dframe
ce_10kb_upstream$strand = '.'
ce_10kb_upstream$score = 1000
ce_10kb_upstream = ce_10kb_upstream[,c('chrom','start','end','geneid','score','strand')]
write.table(ce_10kb_upstream, file = '~/Desktop/ce_10kb_upstream.bed', sep='\t', quote = F, col.names=F, row.names=F)

################################################################################
