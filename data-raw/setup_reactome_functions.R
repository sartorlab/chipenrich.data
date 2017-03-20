#' Build reactome genesets
#'
#' Use the \code{data-raw/NCBI2Reactome_All_Levels.txt} file from \url{http://www.reactome.org/download/current/NCBI2Reactome_All_Levels.txt} to construct the Reactome genesets for the specific \code{org_code}.
#'
#' @param org_code A \code{character} indicating the organism, one of 'Dm', 'Hs', 'Mm', or 'Rn'.
#' @param min_geneset_size An \code{integer} indicating the minimum size for a gene set
#'
#' @return A \code{GeneSet} object representing the Reactome pathways for the particular \code{org_code}.
build_reactome_genesets = function(org_code = c('Dm','Hs','Mm','Rn'), min_geneset_size = 10) {
    org_code = match.arg(org_code)

### {Code bounded needs to be changed per geneset type
    ### Geneset specifics to allow code to generalize
    raw_data_file = 'data-raw/NCBI2Reactome_All_Levels.txt'
    raw_data_columns = c('gene_id','reactome_id','path_url','reactome_name','evidence_code','species')
    raw_data_url = 'http://www.reactome.org/download/current/NCBI2Reactome_All_Levels.txt'
    geneset_file_code = 'reactome'
    doc_geneset_code = 'Reactome'
    doc_geneset_url = 'http://www.reactome.org'
    doc_geneset_accession = 'R-HSA-109688'

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
### Code bounded needs to be changed per geneset type}

    message(sprintf('On %s for %s...', doc_geneset_code, org_name))

    # Read reactome data
    raw_data = read.table(
        file = raw_data_file,
        sep='\t', header=F, quote = '', comment.char = '',
        col.names = raw_data_columns,
        stringsAsFactors=F)

### {Code bounded needs to be changed per geneset type
    # Subset raw_data by species
    species_data = subset(raw_data, species == org_name)

    ### Construct genesets
    # Split the gene_id column by the reactome_id
    genesets = split(species_data$gene_id, species_data$reactome_id)
    # Enforce uniqueness of Entrez IDs
    genesets = lapply(genesets, unique)
    # Enforce min_geneset_size
    genesets = genesets[sapply(genesets,length) >= min_geneset_size]

    ### Construct geneset ID to name mapping
    # Split the reactome_name column by the reactome_id
    geneset_names = split(species_data$reactome_name, species_data$reactome_id)
    # Make unique
    geneset_names = lapply(geneset_names, unique)

    ### Construct all genes
    all_genes = sort(unique(unlist(genesets, use.names = FALSE)))
### Code bounded needs to be changed per geneset type}

    ### Construct the GeneSet object
    gs_obj = new('GeneSet')
    gs_obj@type = doc_geneset_code
    gs_obj@dburl = doc_geneset_url
    gs_obj@organism = org_name
    gs_obj@set.gene = as.environment(genesets)
    gs_obj@all.genes = all_genes
    gs_obj@set.name = as.environment(geneset_names)

    message(sprintf('Saving %s GeneSet object for %s ...', geneset_file_code, org_name))
    eval(parse(text = sprintf("geneset.%s.%s = gs_obj", geneset_file_code, org)))
    save(list = c(sprintf("geneset.%s.%s", geneset_file_code, org)), file = sprintf("data/geneset.%s.%s.RData", geneset_file_code, org), compress='xz')

    message(sprintf('Documenting %s GeneSet object for %s ...', geneset_file_code, org_name))
    # The actual documentation to write, as a roxygen2 block
    # See: http://r-pkgs.had.co.nz/data.html section on Documenting datasets
    doc = c(
        sprintf("#' geneset.%s.%s genesets for %s", geneset_file_code, org, org_name),
        "#'",
        sprintf("#' %s genesets for %s. All genesets are required to have >= %s Entrez IDs.", doc_geneset_code, org_name, min_geneset_size),
        sprintf("#' Built on %s.", date()),
        "#'",
        "#' @format A \\code{GeneSet} object with the following slots:",
        "#' \\describe{",
        sprintf("#'     \\item{type}{A \\code{character} indicating the type of genesets, e.g. %s.}", doc_geneset_code),
        "#'     \\item{dburl}{A \\code{character} of the URL of the database underlying the genesets.}",
        "#'     \\item{organism}{A \\code{character} of the organism, e.g. Homo sapiens.}",
        sprintf("#'     \\item{set.gene}{An \\code{environment} containing a \\code{list} whose keys are database specific accessions (e.g. %s), and whose elements are \\code{character} vectors of Entrez Gene IDs.}", doc_geneset_accession),
        "#'     \\item{all.genes}{A \\code{character} vector of all the Entrez Gene IDs contained over all the genesets in this \\code{type}.}",
        "#'     \\item{set.name}{An \\code{environment} containing a \\code{list} whose keys are database specific accessions, and whose elements are human readable geneset names.}",
        "#' }",
        sprintf("#' @source %s downloaded on 2017-03-19", raw_data_url),
        sprintf('"geneset.%s.%s"', geneset_file_code, org),
        ''
    )

    write(doc, file = 'R/data.R', append = TRUE)
    return(NULL)
}
