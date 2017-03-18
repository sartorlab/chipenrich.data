get_package_version = function(package_name) {
	packages = data.frame(installed.packages()[,c('Package','Version')], row.names = NULL, stringsAsFactors=F)
	version_str = paste(subset(packages, Package == package_name), collapse='_')
	return(version_str)
}

build_GO_genesets = function(org_code = c('Dm','Hs','Mm','Rn'), min_geneset_size = 10) {
	org_code = match.arg(org_code)

	if(org_code == 'Dm') {
		org = 'dme'
	} else if (org_code == 'Hs') {
		org = 'hsa'
	} else if (org_code == 'Mm') {
		org = 'mmu'
	} else if (org_code == 'Rn') {
		org = 'rno'
	}

	message(sprintf('On %s...', org_code))

	require(sprintf('org.%s.eg.db', org_code), character.only = TRUE)
	require(GO.db)

	# Get package versions
	orgdb_version = get_package_version(sprintf('org.%s.eg.db', org_code))
	godb_version = get_package_version('GO.db')

	# Get GO term to Entrez ID mappings from org.*.eg.db
	message('Get GO term to Entrez ID mappings from org.*.eg.db...')
	org_GO_code = sprintf('org.%s.egGO2ALLEGS', org_code)
	org_GO = as.list(get(org_GO_code))

	# Name of organism
	GO_org_name = get(sprintf('org.%s.egORGANISM', org_code))

	# Full names of GO branches
	GO_branch_names = list(
		GOBP = 'Gene Ontology Biological Process',
		GOCC = 'Gene Ontology Cellular Component',
		GOMF = 'Gene Ontology Molecular Function'
	)

	GO_urls = list(
		GOBP = 'http://www.geneontology.org/',
		GOCC = 'http://www.geneontology.org/',
		GOMF = 'http://www.geneontology.org/'
	)

	# Get GO terms for each branch from GO.db
	message('Get GO terms for each branch from GO.db...')
	GO_by_branch = list(
		GOBP = unique(AnnotationDbi::keys(GO.db::GOBPOFFSPRING)),
		GOCC = unique(AnnotationDbi::keys(GO.db::GOCCOFFSPRING)),
		GOMF = unique(AnnotationDbi::keys(GO.db::GOMFOFFSPRING))
	)

	# Filter the mappings from org.*.eg.db by branches from GO.db
	message('Filter the mappings from org.*.eg.db by branches from GO.db...')
	GO_genesets = lapply(GO_by_branch, function(branch_terms){
		org_GO[names(org_GO) %in% branch_terms]
	})

	# Filter by minimum gene set size
	message('Require genesets to have >= 10 genes...')
	GO_genesets = lapply(GO_genesets, function(branch) {
		branch[sapply(branch, length) >= min_geneset_size]
	})

	# Set of all genes in each GO branch
	message('Set of all genes in each GO branch...')
	GO_allgenes = lapply(GO_genesets, function(genesets){
		as.character(sort(Reduce(function(x,y){union(x,y)}, genesets)))
	})

	# Descriptive names of GO terms from GO.db
	message('Descriptive names of GO terms from GO.db...')
	goterms = AnnotationDbi::as.list(GO.db::GOTERM)
	go2desc = Map(function(x) x@Term, goterms)
	GO_names = lapply(GO_by_branch, function(branch_terms){
		go2desc[names(go2desc) %in% branch_terms]
	})

	# Construct a GeneSet object for each GO branch
	construct = lapply(names(GO_names), function(branch){
		message(sprintf('Building %s GeneSet object...', branch))

		gs_obj = new('GeneSet')
		gs_obj@type = GO_branch_names[[branch]]
		gs_obj@dburl = GO_urls[[branch]]
		gs_obj@organism = GO_org_name
		gs_obj@set.gene = as.environment(GO_genesets[[branch]])
		gs_obj@all.genes = GO_allgenes[[branch]]
		gs_obj@set.name = as.environment(GO_names[[branch]])

		message(sprintf('Saving %s GeneSet object...', branch))
		eval(parse(text = sprintf("geneset.%s.%s = gs_obj", branch, org)))
		save(list = c(sprintf("geneset.%s.%s", branch, org)), file = sprintf("data/geneset.%s.%s.RData", branch, org), compress='xz')
		return(NULL)
	})

	for(name in names(GO_branch_names)) {
		message(sprintf('Documenting %s GeneSet object...', name))
		# The actual documentation to write, as a roxygen2 block
		# See: http://r-pkgs.had.co.nz/data.html section on Documenting datasets
		doc = c(
			sprintf("#' geneset.%s.%s genesets for %s", name, org, GO_org_name),
			"#'",
			sprintf("#' %s (%s) genesets for %s. All genesets are required to >= %s Entrez IDs.", GO_branch_names[[name]], name, GO_org_name, min_geneset_size),
			sprintf("#' Built on %s.", date()),
			"#'",
			"#' @format A \\code{GeneSet} object with the following slots:",
			"#' \\describe{",
			"#'     \\item{type}{A \\code{character} indicating the type of genesets, e.g. GOBP.}",
			"#'     \\item{dburl}{A \\code{character} of the URL of the database underlying the genesets.}",
			"#'     \\item{organism}{A \\code{character} of the organism, e.g. Homo sapiens.}",
			"#'     \\item{set.gene}{An \\code{environment} containing a \\code{list} whose keys are database specific accessions (e.g. GO IDs for GO terms), and whose elements are \\code{character} vectors of Entrez Gene IDs.}",
			"#'     \\item{all.genes}{A \\code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \\code{type}.}",
			"#'     \\item{set.name}{An \\code{environment} containing a \\code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}",
			"#' }",
			sprintf("#' @source %s and %s", orgdb_version, godb_version),
			sprintf('"geneset.%s.%s"', name, org),
			''
		)

		write(doc, file = 'R/data.R', append = TRUE)
	}
}
