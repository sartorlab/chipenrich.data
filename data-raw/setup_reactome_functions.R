#' Build reactome genesets
#'
#' Use the \code{data-raw/NCBI2Reactome_All_Levels.txt} file from \url{http://www.reactome.org/download/current/NCBI2Reactome_All_Levels.txt} to construct the Reactome genesets for the specific \code{species}.
#'
#' @param org_code A \code{character} indicating the organism, one of 'Dm', 'Hs', 'Mm', or 'Rn'.
#'
#' @return A \code{GeneSet} object representing the Reactome pathways for the particular \code{org_code}.
build_reactome_genesets = function(org_code = c('Dm','Hs','Mm','Rn'), min_geneset_size = 10) {
    org_code = match.arg(org_code)

    if(org_code == 'Dm') {
        org = 'dme'
        org_name = 'Drosophila melanogaster'
    } else if (org_code == 'Hs') {
        org = 'hsa'
        org_name = 'Homo sapiens'
    } else if (org_code == 'Mm') {
        org = 'mmu'
        org_name = 'Mus musculus'
    } else if (org_code == 'Rn') {
        org = 'rno'
        org_name = 'Rattus norvegicus'
    }

    message(sprintf('On %s...', org_name))

    # Read reactome data
    reactome_raw = read.table(
        file = 'data-raw/NCBI2Reactome_All_Levels.txt',
        sep='\t', header=F, quote = '', comment.char = '',
        col.names = c('gene_id','reactome_id','path_url','reactome_name','evidence_code','species'),
        stringsAsFactors=F)

    # Subset reactome_raw by species
    reactome_species = subset(reactome_raw, species == org_name)

    ### Construct genesets
    # Split the gene_id column by the reactome_id
    genesets = split(reactome_species$gene_id, reactome_species$reactome_id)
    # Enforce uniqueness of Entrez IDs
    genesets = lapply(genesets, unique)
    # Enforce min_geneset_size
    genesets = genesets[sapply(genesets,length) >= min_geneset_size]

    ### Construct geneset ID to name mapping
    # Split the reactome_name column by the reactome_id
    geneset_names = split(reactome_species$reactome_name, reactome_species$reactome_id)
    # Make unique
    geneset_names = lapply(geneset_names, unique)

    ### Construct all genes
    all_genes = sort(unique(unlist(genesets, use.names = FALSE)))

    ### Construct the GeneSet object
    gs_obj = new('GeneSet')
    gs_obj@type = 'Reactome'
    gs_obj@dburl = 'http://www.reactome.org'
    gs_obj@organism = org_name
    gs_obj@set.gene = as.environment(genesets)
    gs_obj@all.genes = all_genes
    gs_obj@set.name = as.environment(geneset_names)

    message('Saving Reactome GeneSet object...')
    eval(parse(text = sprintf("geneset.reactome.%s = gs_obj", org)))
    save(list = c(sprintf("geneset.reactome.%s", org)), file = sprintf("data/geneset.reactome.%s.RData", org), compress='xz')

    message(sprintf('Documenting Reactome GeneSet object for %s ...', org_name))
    # The actual documentation to write, as a roxygen2 block
    # See: http://r-pkgs.had.co.nz/data.html section on Documenting datasets
    doc = c(
        sprintf("#' geneset.reactome.%s genesets for %s", org, org_name),
        "#'",
        sprintf("#' Reactome genesets for %s. All genesets are required to have >= %s Entrez IDs.", org_name, min_geneset_size),
        sprintf("#' Built on %s.", date()),
        "#'",
        "#' @format A \\code{GeneSet} object with the following slots:",
        "#' \\describe{",
        "#'     \\item{type}{A \\code{character} indicating the type of genesets, e.g. Reactome.}",
        "#'     \\item{dburl}{A \\code{character} of the URL of the database underlying the genesets.}",
        "#'     \\item{organism}{A \\code{character} of the organism, e.g. Homo sapiens.}",
        "#'     \\item{set.gene}{An \\code{environment} containing a \\code{list} whose keys are database specific accessions (e.g. R-HSA-109688), and whose elements are \\code{character} vectors of Entrez Gene IDs.}",
        "#'     \\item{all.genes}{A \\code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \\code{type}.}",
        "#'     \\item{set.name}{An \\code{environment} containing a \\code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}",
        "#' }",
        "#' @source http://www.reactome.org/download/current/NCBI2Reactome_All_Levels.txt downloaded on 2017-03-19",
        sprintf('"geneset.reactome.%s"', org),
        ''
    )

    write(doc, file = 'R/data.R', append = TRUE)
    return(NULL)
}
