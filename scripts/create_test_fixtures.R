suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
})

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

sce <- SingleCellExperiment(
  assays = list(
    counts = counts,
    logcounts = log1p(counts),
    dense_logcounts = as.matrix(log1p(counts))
  ),
  colData = DataFrame(
    cell_type = factor(c("T cell", "B cell", "T cell", "B cell")),
    batch = c("one", "one", "two", "two"),
    score = c(0.1, 0.2, 0.3, 0.4),
    row.names = colnames(counts)
  )
)

saveRDS(sce, "testdata/sce.rds", compress = "gzip")

groups <- colData(sce)$cell_type
expected <- sapply(levels(groups), function(group) {
  log1p(rowMeans(expm1(assay(sce, "logcounts")[, groups == group, drop = FALSE])))
})
write.table(
  cbind(gene = rownames(sce), as.data.frame(expected, check.names = FALSE)),
  "testdata/sce.expected.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
