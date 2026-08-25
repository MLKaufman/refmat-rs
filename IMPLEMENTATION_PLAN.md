# Refmat implementation plan

## Objective

Build a standalone Rust CLI that automatically reads the three common
single-cell containers and produces the same genes-by-groups reference matrix:

```text
refmat head <seurat.rds|sce.rds|data.h5ad>
refmat build <input> --column <annotation>
```

R, Bioconductor, Python, and Scanpy/AnnData are development-only dependencies.
The installed executable must not require any of them.

## Format detection and defaults

Detection is content-based:

1. The HDF5 signature identifies H5AD, after which the AnnData root encoding and
   required `obs`/`var` groups are validated.
2. Other inputs are parsed as RDS and classified from the root S4 class as
   `Seurat` or `SingleCellExperiment`.

Format-aware defaults keep the short command useful:

| Input | Metadata | Default matrix | Orientation |
| --- | --- | --- | --- |
| Seurat v3/v4 RDS | `meta.data` | active assay, `data` slot | features × cells |
| Seurat v5 RDS | `meta.data` | active assay, `data` layer | features × cells |
| SingleCellExperiment RDS | `colData` | `logcounts`, then `counts`, then first assay | features × cells |
| AnnData H5AD | `obs` | `X` | cells × features |

`--assay` overrides an R assay and `--layer` overrides a Seurat or AnnData
layer. `--scale auto|log1p|linear` controls how stored values are interpreted.
Automatic scale selection treats a matrix named `counts` as linear and all
other names as log1p-normalized. AnnData itself does not specify the semantic
scale of `X`, so users must select `--scale linear` when `X` contains counts.

## Shared calculation

Every adapter supplies ordered cell metadata, ordered feature identifiers,
matrix dimensions and orientation, and a bounded-memory stream of dense values
or sparse entries.

For each group, refmat calculates:

```text
log1p(mean(linear expression))
```

Log1p input is converted with `expm1` before accumulation. Implicit sparse
zeroes remain in the group-size divisor. Missing annotations are excluded.
Accumulation and output use `f64`.

## RDS adapters

The RDS reader keeps large vectors lazy and decompresses compressed RDS data to
a temporary backing file. Only selected metadata and matrix ranges are read.

### Seurat

The adapter reads `meta.data`, `active.assay`, `assays`, and the stored object
version. It detects each selected assay structurally. Legacy Seurat v3/v4
`Assay` objects expose `counts`, `data`, and `scale.data` slots with matrix
dimnames; Seurat v5 `Assay5` objects expose named layers and `LogMap` cell and
feature membership maps. Legacy slots accept `dgCMatrix` or dense numeric
matrices, while Assay5 layers currently accept `dgCMatrix`. Both layouts require
a full matrix in metadata cell order.

### SingleCellExperiment

The adapter reads the standard inherited `colData` and `assays` slots. It
supports ordinary in-memory Matrix `dgCMatrix` assays and dense integer/double R
matrices. It deliberately rejects `DelayedArray`, `HDF5Array`, and other
external-backed representations because those serialized objects refer to
additional R classes, files, or native pointers.

## H5AD adapter

HDF5 is statically built into the Rust binary with gzip and LZF filter support;
the RDS xz decoder's liblzma is also linked statically.
The adapter directly implements the stable AnnData on-disk encodings needed by
this workflow:

- `obs` and `var` dataframe indices and column ordering;
- categorical, string, numeric, and nullable numeric metadata columns;
- dense numeric arrays;
- CSR and CSC sparse groups containing `data`, `indices`, and `indptr`.

AnnData matrix orientation is transposed logically during aggregation; the
matrix is not physically transposed. Dense data is read by row blocks. Sparse
indices and values are read in one-million-entry chunks, while `indptr` remains
the only complete sparse vector held in memory.

## CLI contract

```text
refmat inspect <FILE> [--depth N] [--full]
refmat head <FILE> [-n ROWS]
refmat col <FILE> <COLUMN>
refmat build <FILE> --column NAME [--assay NAME] [--layer NAME]
             [--scale auto|log1p|linear] [--output FILE]
```

The default output is `<input-stem>.refmat.tsv`. The input filename is always
required; refmat never selects an arbitrary object from the working directory.

## Validation completed

- The supplied 1.9 GB Seurat v5 object generated a 32,285 × 18 result matching
  the equivalent R/Matrix calculation to a maximum absolute error of
  `8.881784e-16`.
- Structurally faithful Seurat v3 and v4 fixtures validate stored-version and
  assay-layout detection, metadata display, grouping counts, sparse `counts` and
  `data`, and dense `scale.data`.
- A genuine Bioconductor `SingleCellExperiment` fixture validates `colData`,
  factor ordering, sparse `dgCMatrix`, dense assays, `logcounts`, and `counts`.
- A genuine Python AnnData fixture validates categorical/string/numeric `obs`,
  dense log1p `X`, CSR and CSC sparse matrices, gzip and LZF compression, and a
  `counts` layer.
- Rust unit tests, formatting, and Clippy run without warnings.

## Remaining hardening

1. Add fixtures for AnnData nullable columns, missing categories, and float32
   data.
2. Add SCE fixtures with absent dimnames, integer dense counts, and additional
   Bioconductor release versions.
3. Align cells by identifier rather than requiring positional equality.
4. Add Seurat split-layer and partial-layer support.
5. Add output allocation limits, malformed-file tests, and atomic output writes.
6. Publish prebuilt macOS/Linux binaries and checksums so end users do not need
   the Rust/CMake build toolchain.

## Distribution

`cargo install --git https://github.com/MLKaufman/refmat-rs.git --locked` builds a
static-HDF5 executable from the repository. Building requires Rust 1.85+, CMake,
and a C compiler, but running the resulting executable has no R, Python, or
system-HDF5 dependency.
