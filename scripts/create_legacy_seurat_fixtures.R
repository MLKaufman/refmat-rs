suppressPackageStartupMessages(library(Matrix))

# These minimal S4 definitions reproduce the serialized slots used by Seurat
# v3/v4 without requiring old Seurat releases in the development environment.
setClass(
  "Assay",
  slots = c(
    counts = "ANY",
    data = "ANY",
    scale.data = "ANY",
    assay.orig = "character",
    var.features = "vector",
    meta.features = "data.frame",
    misc = "list",
    key = "character"
  )
)
setClass(
  "Seurat",
  slots = c(
    assays = "list",
    meta.data = "data.frame",
    active.assay = "character",
    active.ident = "factor",
    graphs = "list",
    neighbors = "list",
    reductions = "list",
    images = "list",
    project.name = "character",
    misc = "list",
    version = "ANY",
    commands = "list",
    tools = "list"
  )
)

counts <- Matrix(
  c(
    1, 0, 5, 0,
    0, 2, 0, 4,
    3, 0, 1, 0
  ),
  nrow = 3,
  byrow = TRUE,
  sparse = TRUE,
  dimnames = list(c("GeneA", "GeneB", "GeneC"), paste0("Cell", 1:4))
)
normalized <- log1p(counts)
metadata <- data.frame(
  cell_type = factor(c("T cell", "B cell", "T cell", "B cell")),
  batch = c("one", "one", "two", "two"),
  score = c(0.1, 0.2, 0.3, 0.4),
  row.names = colnames(counts)
)
assay <- new(
  "Assay",
  counts = counts,
  data = normalized,
  scale.data = as.matrix(normalized),
  assay.orig = character(),
  var.features = character(),
  meta.features = data.frame(row.names = rownames(counts)),
  misc = list(),
  key = "rna_"
)

make_fixture <- function(version, path) {
  object <- new(
    "Seurat",
    assays = list(RNA = assay),
    meta.data = metadata,
    active.assay = "RNA",
    active.ident = metadata$cell_type,
    graphs = list(),
    neighbors = list(),
    reductions = list(),
    images = list(),
    project.name = "refmat fixture",
    misc = list(),
    version = package_version(version),
    commands = list(),
    tools = list()
  )
  saveRDS(object, path, compress = "gzip")
}

make_fixture("3.2.3", "testdata/seurat-v3.rds")
make_fixture("4.4.0", "testdata/seurat-v4.rds")

groups <- metadata$cell_type
expected <- sapply(levels(groups), function(group) {
  log1p(rowMeans(expm1(normalized[, groups == group, drop = FALSE])))
})
write.table(
  cbind(gene = rownames(counts), as.data.frame(expected, check.names = FALSE)),
  "testdata/seurat-legacy.expected.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
