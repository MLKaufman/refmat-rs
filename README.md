# refmat

`refmat` is an experimental standalone Rust CLI for inspecting single-cell
objects and generating genes-by-cell-type reference matrices without requiring
R or Python at runtime.

Input format and object type are detected automatically. The current prototype
supports:

- Seurat v5 objects in RDS files with in-memory `Assay5` `dgCMatrix` layers
- Bioconductor `SingleCellExperiment` objects in RDS files with in-memory dense
  numeric or `dgCMatrix` assays
- Scanpy/AnnData H5AD files with dense, CSR, or CSC `X`/layers, including gzip
  and LZF compression

## Build

```bash
cargo build --release
```

The executable is `target/release/refmat`.

Building requires Rust 1.85 or newer, CMake, and a C compiler. HDF5, gzip/LZF,
and the RDS xz decoder are built into the executable; system HDF5 and liblzma
installations are not required.

## Install

From a local clone of this repository:

```bash
cargo install --path .
```

This builds an optimized binary and installs `refmat` into Cargo's binary
directory, normally `~/.cargo/bin`. Make sure that directory is included in
your `PATH`, then confirm the installation:

```bash
refmat --version
```

Install directly from GitHub without cloning the repository first:

```bash
cargo install --git https://github.com/MLKaufman/refmat.git --locked
```

The `--locked` option uses the dependency versions recorded in `Cargo.lock` for
a reproducible build. Re-run the same command with `--force` to replace an
existing installation:

```bash
cargo install --git https://github.com/MLKaufman/refmat.git --locked --force
```

## Commands

Inspect an object without loading its expression matrix into memory:

```bash
refmat inspect testdata/so.rds
refmat inspect sample.h5ad
```

Print cell metadata as a readable, aligned table:

```bash
refmat head testdata/so.rds
refmat head testdata/so.rds -n 20
refmat head sce.rds
refmat head sample.h5ad
```

Count cells for each value in a cell metadata column:

```bash
refmat col testdata/so.rds celltypes
refmat col sce.rds cell_type
refmat col sample.h5ad cell_type
```

Build a reference matrix from the active assay's normalized `data` layer:

```bash
refmat build testdata/so.rds --column celltypes
```

The same command works for a `SingleCellExperiment` or AnnData file:

```bash
refmat build sce.rds --column cell_type
refmat build sample.h5ad --column cell_type
```

Defaults are format-aware:

- Seurat: active assay and its `data` layer
- SingleCellExperiment: `logcounts`, then `counts`, then the first assay
- AnnData: `X`

The default output is `testdata/so.refmat.tsv`. Assay, layer, and output can be
selected explicitly:

```bash
refmat build testdata/so.rds \
  --column celltypes \
  --assay RNA \
  --layer data \
  --output testdata/custom.refmat.tsv
```

Select a SingleCellExperiment assay or AnnData layer explicitly:

```bash
refmat build sce.rds --column cell_type --assay logcounts
refmat build sample.h5ad --column cell_type --layer counts
```

For `data`, values are calculated as:

```text
log1p(mean(expm1(log_normalized_expression)))
```

For `counts`, values are calculated as `log1p(mean(counts))`.

AnnData does not record whether `X` is raw or log-normalized. Automatic mode
assumes `X` and non-`counts` assay/layer names are log1p-normalized. Override
that interpretation when necessary:

```bash
refmat build raw-counts.h5ad --column cell_type --scale linear
refmat build sample.h5ad --column cell_type --layer normalized --scale log1p
```

The file contents, not the filename extension, determine whether an input is
H5AD or RDS. RDS inputs are then classified from their S4 class.

## Validation on the supplied fixture

The supplied `testdata/so.rds` contains:

- 61,844 cells
- 32,285 features
- `counts`, `data`, and `scale.data` layers
- 59 cell metadata columns
- an 18-level `celltypes` annotation
- 135,265,589 nonzero entries in the `data` layer

The release build generated the complete 32,285 × 18 reference matrix in about
16 seconds on the development machine, including gzip decompression and RDS
structure parsing.

Comparison with the equivalent R/Matrix calculation produced:

```text
max absolute error:  8.881784e-16
mean absolute error: 2.137444e-18
all.equal tolerance 1e-12: TRUE
```

Re-run that comparison with:

```bash
Rscript scripts/validate_reference.R \
  testdata/so.rds \
  testdata/so.refmat.tsv \
  celltypes
```

R is used only for development validation and is not linked into the binary.

## Current limitations

- Seurat expression layers must currently be `dgCMatrix`; SCE supports dense
  numeric matrices and `dgCMatrix`.
- Metadata-to-layer alignment currently requires the layer to contain all cells
  in Seurat metadata order.
- Split layers are not joined automatically.
- BPCells, DelayedArray/HDF5Array inside RDS, and external-pointer matrix
  backends are not supported. Export those objects to an ordinary in-memory
  assay or H5AD first.
- H5AD `X` and named layers are supported, but `.raw/X` is not yet selectable.
- H5AD annotation columns used for grouping must be categorical or string.
- Character annotations preserve first-seen order; factor annotations preserve
  R factor-level order; AnnData categoricals preserve category order.
- Cells with missing annotations are excluded with no output column.
- Each invocation decompresses a gzip RDS into a temporary backing file. A
  persistent cache could improve repeated-command latency.

See [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) for the architecture,
milestones, test matrix, and distribution plan.
