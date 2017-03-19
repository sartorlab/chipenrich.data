source('data-raw/setup_reactome_functions.R')
devtools::load_all()

build_reactome_genesets(org_code = 'Hs', min_geneset_size = 10)
build_reactome_genesets(org_code = 'Mm', min_geneset_size = 10)
build_reactome_genesets(org_code = 'Rn', min_geneset_size = 10)
build_reactome_genesets(org_code = 'Dm', min_geneset_size = 10)
