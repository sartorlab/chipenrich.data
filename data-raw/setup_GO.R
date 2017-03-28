source('data-raw/setup_GO_functions.R')
devtools::load_all()

build_GO_genesets(org_code = 'Hs', min_geneset_size = 10)
build_GO_genesets(org_code = 'Mm', min_geneset_size = 10)
build_GO_genesets(org_code = 'Rn', min_geneset_size = 10)
build_GO_genesets(org_code = 'Dm', min_geneset_size = 10)
build_GO_genesets(org_code = 'Dr', min_geneset_size = 10)
