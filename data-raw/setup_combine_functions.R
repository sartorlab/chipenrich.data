build_combine = function(type = c('BioCarta', 'BioGrid', 'CTD', 'DrugBank', 'Hallmark', 'Immunologic', 'KEGG', 'MicroRNA', 'Oncogenic', 'pfam', 'TF'), org_code = c('Hs'), min_geneset_size = 10) {
    org_code = match.arg(org_code)

    if(org_code == 'Dm') {
        org = 'dme'
        organism = 'Drosophila melanogaster'
    } else if (org_code == 'Hs') {
        org = 'hsa'
        organism = 'Homo sapiens'
    } else if (org_code == 'Mm') {
        org = 'mmu'
        organism = 'Mus musculus'
    } else if (org_code == 'Rn') {
        org = 'rno'
        organism = 'Rattus norvegicus'
    } else if (org_code == 'Dr') {
        org = 'dre'
        organism = 'Danio rerio'
    }

    type = match.arg(type)

    if(type == 'BioCarta') {
        file = 'data-raw/BioCarta_combine_file.txt'
        geneset_name = 'BioCarta'
        geneset_code = 'biocarta_pathway'
        db_url = 'https://cgap.nci.nih.gov/Pathways/BioCarta_Pathways'
    } else if (type == 'CTD') {
        file = 'data-raw/CTD_combine_file_hsa.txt'
        geneset_name = 'Comparative Toxicogenomics Database'
        geneset_code = 'ctd'
        db_url = 'http://ctdbase.org'
    } else if (type == 'DrugBank') {
        file = 'data-raw/DrugBank_combine_file.txt'
        geneset_name = 'DrugBank'
        geneset_code = 'drug_bank'
        db_url = 'https://www.drugbank.ca'
    } else if (type == 'Hallmark') {
        file = 'data-raw/Hallmark_combine_file.txt'
        geneset_name = 'Hallmark (MSigDB)'
        geneset_code = 'hallmark'
        db_url = 'http://software.broadinstitute.org/gsea/msigdb/collections.jsp#H'
    } else if (type == 'Immunologic') {
        file = 'data-raw/Immunologic_Signatures_combine_file.txt'
        geneset_name = 'Immunologic Signatures (MSigDB)'
        geneset_code = 'immunologic'
        db_url = 'http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C7'
    } else if (type == 'KEGG') {
        file = 'data-raw/KEGG_combine_file.txt'
        geneset_name = 'KEGG Pathways'
        geneset_code = 'kegg_pathway'
        db_url = 'http://kegg.jp'
    } else if (type == 'MicroRNA') {
        file = 'data-raw/MicroRNA_targets_combine_file.txt'
        geneset_name = 'MicroRNA Targets (MSigDB)'
        geneset_code = 'microrna'
        db_url = 'http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C3'
    } else if (type == 'Oncogenic') {
        file = 'data-raw/Oncogenic_Signatures_combine_file.txt'
        geneset_name = 'Oncogenic Signatures (MSigDB)'
        geneset_code = 'oncogenic'
        db_url = 'http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C6'
    } else if (type == 'pfam') {
        file = 'data-raw/pfam_combine_file.txt'
        geneset_name = 'Pfam'
        geneset_code = 'pfam'
        db_url = 'http://pfam.xfam.org'
    } else if (type == 'TF') {
        file = 'data-raw/Transcription_factor_targets_combine_file.txt'
        geneset_name = 'Transcription Factor Targets (MSigDB)'
        geneset_code = 'transcription_factors'
        db_url = 'http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C3'
    } else if (type == 'BioGrid') {
        file = 'data-raw/BIOgrid_combine_file_hsa.txt'
        geneset_name = 'BioGRID Protein Interactions'
        geneset_code = 'protein_interaction_biogrid'
        db_url = 'https://thebiogrid.org'
    }

    message(sprintf('On %s...', org_code))

    ### Read data
    raw_data = read.table(file, header = T, sep='\t', quote = '\"', as.is = TRUE)

    ### Build the constituent parts
    # Want these two remain as character, do not want Entrez IDs to be integers or numerics
    genesets = lapply(raw_data$EntrezID, function(genes){
        unlist(strsplit(genes, ','))
    })
    names(genesets) = raw_data$GenesetID

    allgenes = as.character(unlist(genesets))

    gs_names = as.list(raw_data$GenesetName)
    names(gs_names) = raw_data$GenesetID

    ### Build the GeneSet object
    message(sprintf('Building %s GeneSet object...', geneset_code))

    gs_obj = new('GeneSet')
    gs_obj@type = geneset_name
    gs_obj@dburl = db_url
    gs_obj@organism = organism
    gs_obj@set.gene = as.environment(genesets)
    gs_obj@all.genes = allgenes
    gs_obj@set.name = as.environment(gs_names)

    ### Save it
    message(sprintf('Saving %s GeneSet object...', geneset_code))
    eval(parse(text = sprintf("geneset.%s.%s = gs_obj", geneset_code, org)))
    save(list = c(sprintf("geneset.%s.%s", geneset_code, org)), file = sprintf("data/geneset.%s.%s.RData", geneset_code, org), compress='xz')

    ### Document it
    message(sprintf('Documenting %s GeneSet object...', geneset_name))
    # The actual documentation to write, as a roxygen2 block
    # See: http://r-pkgs.had.co.nz/data.html section on Documenting datasets
    doc = c(
        sprintf("#' geneset.%s.%s genesets for %s", geneset_code, org, geneset_name),
        "#'",
        sprintf("#' %s (%s) genesets. All genesets are required to have >= %s Entrez IDs.", geneset_name, geneset_code, min_geneset_size),
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
        sprintf("#' @source %s", db_url),
        sprintf('"geneset.%s.%s"', geneset_code, org),
        ''
    )

    write(doc, file = 'R/genesets_combine_doc.R', append = TRUE)
    return(NULL)
}
