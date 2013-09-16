# S4 class for storing genesets. 
setClass("GeneSet",representation(
  set.gene = "environment",
  type = "character",
  set.name = "environment",
  all.genes = "character",
  organism = "character",
  dburl = "character"
),prototype(
  set.gene = new.env(parent=emptyenv()),
  type = "",
  set.name = new.env(parent=emptyenv()),
  all.genes = "",
  organism = "",
  dburl = ""
),
  package = "chipenrich.data"
);

# S4 class for storing locus definitions. 
setClass("LocusDefinition",representation(
  dframe = "data.frame",
  granges = "GenomicRanges",
  chrom2iranges = "list",
  genome.build = "character",
  organism = "character"
),
  package = "chipenrich.data"
);
