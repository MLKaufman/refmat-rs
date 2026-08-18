args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: Rscript scripts/validate_reference.R <seurat.rds> <refmat.tsv> <column>",
    call. = FALSE
  )
}

input_path <- args[[1L]]
reference_path <- args[[2L]]
column <- args[[3L]]

so <- readRDS(input_path)
assay <- so@assays[[so@active.assay]]
expression <- assay@layers[["data"]]
groups <- so@meta.data[[column]]

if (!is.factor(groups)) {
  groups <- factor(groups, levels = unique(groups))
}

group_names <- levels(groups)
feature_names <- attr(assay@features, "dimnames")[[1L]]
rust <- read.delim(reference_path, row.names = 1L, check.names = FALSE)

stopifnot(
  identical(colnames(rust), group_names),
  identical(rownames(rust), feature_names),
  nrow(rust) == nrow(expression)
)

expected <- vapply(
  group_names,
  function(group) {
    log1p(Matrix::rowMeans(expm1(expression[, groups == group, drop = FALSE])))
  },
  numeric(nrow(expression))
)

difference <- abs(as.matrix(rust) - expected)
cat("dimensions\t", paste(dim(rust), collapse = "x"), "\n", sep = "")
cat("max_abs_error\t", format(max(difference), scientific = TRUE), "\n", sep = "")
cat("mean_abs_error\t", format(mean(difference), scientific = TRUE), "\n", sep = "")
cat(
  "all_equal_1e-12\t",
  isTRUE(all.equal(
    as.matrix(rust),
    expected,
    tolerance = 1e-12,
    check.attributes = FALSE
  )),
  "\n",
  sep = ""
)

