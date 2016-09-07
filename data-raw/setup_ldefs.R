library(GenomicRanges)
devtools::load_all()
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

################################################################################
################################################################################
# Object setup
################################################################################
################################################################################

#######################################
# Get all the transcripts
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
tx_gr = transcripts(txdb, columns = c('GENEID','TXNAME'))
seqinfo(tx_gr) = seqinfo(txdb)

# Get transcripts that have valid GENEIDs (i.e. not NA)
tx_gr = tx_gr[!is.na(as.character(mcols(tx_gr)$GENEID))]

# Destrand tx_gr and sort it
strand(tx_gr) = '*'
tx_gr = sort(tx_gr)

#######################################
# Construct a genome GRanges
n_seqnames = length(seqnames(seqinfo(txdb)))
genome_gr = GenomicRanges::GRanges(
	seqnames = seqnames(seqinfo(txdb)),
	ranges = IRanges::IRanges(start = rep.int(1, n_seqnames), end = seqlengths(seqinfo(txdb))),
)
seqinfo(genome_gr) = seqinfo(txdb)

#######################################
# Get all the TSSs and make them unique (collapsing GENEIDs into IntegerList)
# Treat geneid as integers so that splitAsList works faster
tss_gr = GenomicRanges::GRanges(
	seqnames = seqnames(tx_gr),
	ranges = IRanges::IRanges(start = start(tx_gr), end = start(tx_gr)),
	geneid = as.integer(mcols(tx_gr)$GENEID)
)
seqinfo(tss_gr) = seqinfo(tx_gr)
tss_gr = sort(tss_gr)

# Get the unique TSSs and add the GENEIDs back
utss_gr = unique(granges(tss_gr))
mcols(utss_gr)$geneid = splitAsList(tss_gr$geneid, match(tss_gr, utss_gr))
mcols(utss_gr)$geneid = IntegerList(lapply(mcols(utss_gr)$geneid, unique))

#######################################
# Make tx_gr unique and use it for nearest gene definition

utx_gr = unique(granges(tx_gr))
mcols(utx_gr)$geneid = splitAsList(as.integer(tx_gr$GENEID), match(tx_gr, utx_gr))
mcols(utx_gr)$geneid = IntegerList(lapply(mcols(utx_gr)$geneid, unique))

#######################################
# Get all the exons and introns
exon_gr = exons(txdb, columns = c('GENEID','TXNAME'))
exon_gr = exon_gr[sapply(exon_gr$GENEID, length) > 0]

intron_gr = GenomicFeatures::intronsByTranscript(txdb)

################################################################################
################################################################################
# Determine quantities for creation of nearest_tss and x kb definitions
################################################################################
################################################################################

#######################################
# Get distances of preceding and following TSSs

# Split into GRangesList so that follow and precede stay intra chromosome
# Keep only the chromosomes that have two or more TSSs.
utss_split_gr = split(utss_gr, seqnames(utss_gr))
utss_split_gr = utss_split_gr[sapply(utss_split_gr, length) >= 2]

# Determine precede and follow: idxs, dists, midpoints, 1k, 5k, 10k
# NOTE: Since strand is *, precede means upstream (left),
# and follow means downstream (right)
utss_split_gr = endoapply(utss_split_gr, function(chr){
	# This looks wrong, but is not. I'm confused too.
	chr$precede_idx = GenomicRanges::follow(chr, chr)
	chr$follow_idx = GenomicRanges::precede(chr, chr)

	# Set first and last indexes for precede and follow to first and last
	# indices, respectively so that vectorized operations work
	chr$precede_idx[1] = 1
	chr$follow_idx[length(chr)] = length(chr)

	# Determine distances
	chr$precede_dist = GenomicRanges::distance(
		x = chr,
		y = chr[chr$precede_idx]
	)
	chr$follow_dist = GenomicRanges::distance(
		x = chr,
		y = chr[chr$follow_idx]
	)

	# Reset the preceding and following distances for the first and last transcripts
	# to the distances to the beginning and end of the chromosome, respectively.
	chr$precede_dist[1] = start(chr)[1] - 1
	chr$follow_dist[length(chr)] = end(genome_gr[seqnames(genome_gr) == unique(seqnames(chr))]) - end(chr)[length(chr)]

	# Determine midpoints between preceeding and following TSSs based on adjusted distances
	chr$precede_mid = floor(start(ranges(chr)) - (chr$precede_dist / 2))
	chr$follow_mid = floor(start(ranges(chr)) + (chr$follow_dist / 2))

	# The nearest TSS definition is the easiest as it is the midpoints between TSSs
	chr$start_tss = chr$precede_mid
	chr$end_tss = chr$follow_mid

	# Determine the start and end of the x kb locus definitions
	for(ldf in c(1000,5000,10000)) {
		# Determine ldf upstream and downstream of the TSSs
		upstream = start(ranges(chr)) - ldf
		downstream = start(ranges(chr)) + ldf

		# Legend
		# (TSS) is the TSS of interest
		# mid is the midpoint between adjacent TSSs
		# (start = precede_mid or end = follow_mid)
		# u is the extension x bp upstream from the (TSS)
		# d is the extension x bp downstream from the (TSS)

		# Determine where to start ldf
		# --TSS-------------------mid----------u--------(TSS) (pick u)
		# --TSS--------u----------mid-------------------(TSS) (pick mid)
		start = chr$precede_mid
		replace = which(start < upstream)
		start[replace] = upstream[replace]

		# Determine where to end ldf
		# --(TSS)-------------------mid----------d--------TSS (pick mid)
		# --(TSS)--------d----------mid-------------------TSS (pick d)
		end = chr$follow_mid
		replace = which(end > downstream)
		end[replace] = downstream[replace]

		# Combine into a data.frame to append to mcols(chr)
		tmp_df = data.frame(
			upstream,
			downstream,
			start,
			end,
			stringsAsFactors=F
		)
		colnames(tmp_df) = c(
			sprintf('upstream_%s', ldf),
			sprintf('downstream_%s', ldf),
			sprintf('start_%s', ldf),
			sprintf('end_%s', ldf)
		)

		# Add the columns
		mcols(chr) = cbind(mcols(chr), tmp_df)
	}

	return(chr)
})

#######################################
# Recombine utss_split_gr
utss_gr = Reduce(c, utss_split_gr)

################################################################################
################################################################################
# Build GRanges for locus definitions
################################################################################
################################################################################

ldef_tss_gr = GRanges(
	seqnames = seqnames(utss_gr),
	ranges = IRanges(start = utss_gr$start_tss, end = utss_gr$end_tss),
	geneid = utss_gr$geneid
)
ldef_1kb_gr = GRanges(
	seqnames = seqnames(utss_gr),
	ranges = IRanges(start = utss_gr$start_1000, end = utss_gr$end_1000),
	geneid = utss_gr$geneid
)
ldef_5kb_gr = GRanges(
	seqnames = seqnames(utss_gr),
	ranges = IRanges(start = utss_gr$start_5000, end = utss_gr$end_5000),
	geneid = utss_gr$geneid
)
ldef_10kb_gr = GRanges(
	seqnames = seqnames(utss_gr),
	ranges = IRanges(start = utss_gr$start_10000, end = utss_gr$end_10000),
	geneid = utss_gr$geneid
)

################################################################################
################################################################################
# Build data.frames for locus definitions
################################################################################
################################################################################

# ldef_tss_df =
# ldef_1kb_df =
# ldef_5kb_df =
# ldef_10kb_df =
