## ---- eval=FALSE---------------------------------------------------------
#  methods::setClass("LocusDefinition", methods::representation(
#    dframe = "data.frame",
#    granges = "GRanges",
#    genome.build = "character",
#    organism = "character"
#  ),
#    package = "chipenrich.data"
#  );

## ---- eval=FALSE---------------------------------------------------------
#  # S4 class for storing genesets.
#  methods::setClass("GeneSet",representation(
#    set.gene = "environment",
#    type = "character",
#    set.name = "environment",
#    all.genes = "character",
#    organism = "character",
#    dburl = "character"
#  ), methods::prototype(
#    set.gene = new.env(parent=emptyenv()),
#    type = "",
#    set.name = new.env(parent=emptyenv()),
#    all.genes = "",
#    organism = "",
#    dburl = ""
#  ),
#    package = "chipenrich.data"
#  )

