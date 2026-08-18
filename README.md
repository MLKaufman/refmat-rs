# refmat

`refmat` is an experimental standalone Rust CLI for inspecting Seurat RDS files
and generating genes-by-cell-type reference matrices without requiring R at
runtime.

The current prototype supports the supplied Seurat v5 fixture and ordinary
in-memory `Assay5` `dgCMatrix` layers.

## Build

```bash
cargo build --release
```

The executable is `target/release/refmat`.

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

Once the repository is hosted on GitHub, it can be installed directly without
cloning it first:

```bash
cargo install --git https://github.com/OWNER/refmat.git --locked
```

Replace `OWNER` with the GitHub account or organization hosting the repository.
The `--locked` option uses the dependency versions recorded in `Cargo.lock` for
a reproducible build. Re-run the same command with `--force` to replace an
existing installation:

```bash
cargo install --git https://github.com/OWNER/refmat.git --locked --force
```

## Commands

Inspect the generic Seurat/R object structure without loading large vectors:

```bash
refmat inspect testdata/so.rds
```

Print cell metadata as TSV:

```bash
refmat head testdata/so.rds
refmat head testdata/so.rds -n 20
```

Build a reference matrix from the active assay's normalized `data` layer:

```bash
refmat build testdata/so.rds --column celltypes
```

The default output is `testdata/so.refmat.tsv`. Assay, layer, and output can be
selected explicitly:

```bash
refmat build testdata/so.rds \
  --column celltypes \
  --assay RNA \
  --layer data \
  --output testdata/custom.refmat.tsv
```

For `data`, values are calculated as:

```text
log1p(mean(expm1(log_normalized_expression)))
```

For `counts`, values are calculated as `log1p(mean(counts))`.

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

- Only `dgCMatrix` expression layers are accepted by `build`.
- Metadata-to-layer alignment currently requires the layer to contain all cells
  in Seurat metadata order.
- Split layers are not joined automatically.
- BPCells, HDF5, DelayedArray, and external-pointer matrix backends are not
  supported.
- Character annotations preserve first-seen order; factor annotations preserve
  R factor-level order.
- Cells with missing annotations are excluded with no output column.
- Each invocation decompresses a gzip RDS into a temporary backing file. A
  persistent cache could improve repeated-command latency.

See [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) for the architecture,
milestones, test matrix, and distribution plan.
