# Refmat implementation plan

## Objective

Build `refmat`, a standalone Rust CLI that reads ordinary in-memory Seurat v4/v5
objects stored as RDS files without requiring R at runtime. The first supported
workflow is:

```text
refmat head so.rds
refmat build so.rds --column cell_types
```

`build` writes a dense, tab-separated genes-by-groups reference matrix. Its
default output is `<input-stem>.refmat.tsv`.

R, Seurat, and clustifyr are development-only dependencies used to create test
fixtures and expected outputs.

## Supported initial scope

- RDS serialization readable by the pinned `rds2rust` version.
- Seurat v4 `Assay` objects with `data` and `counts` slots.
- Seurat v5 `Assay5` objects with one explicitly resolvable `data` or `counts`
  layer.
- Dense numeric matrices and Matrix `dgCMatrix` expression layers.
- Character, factor, integer, double, and logical cell metadata columns.
- Mean aggregation of normalized log1p expression, matching the common
  clustifyr `seurat_ref(..., method = "mean", if_log = TRUE)` path.

Initially unsupported inputs must fail with an actionable error. These include
BPCells, HDF5-backed matrices, DelayedArray objects, external pointers, and
ambiguous split layers that need Seurat's `JoinLayers` behavior.

## CLI contract

```text
refmat inspect <FILE>
refmat head <FILE> [-n <ROWS>]
refmat metadata <FILE> --columns
refmat build <FILE> --column <NAME> [--assay <NAME>] [--layer <NAME>]
             [--output <FILE>] [--input-log | --input-linear]
```

Defaults:

- assay: the Seurat object's `active.assay`
- layer: `data`; no silent fallback once the user supplies `--layer`
- aggregation: mean
- input interpretation: log1p-normalized for `data`, linear for `counts`
- output transform: `log1p(mean(linear expression))`
- output: `<input-stem>.refmat.tsv`

The input filename remains mandatory. The program must never choose an arbitrary
RDS file from the current directory.

## Architecture

```text
CLI
  -> RDS reader adapter
       -> generic R object/path navigation
            -> Seurat v4/v5 adapter
                 -> MetadataTable + ExpressionLayer
                      -> grouped averaging
                           -> atomic TSV writer
```

Only the RDS and Seurat adapters may depend on `rds2rust` object details. The
averaging code operates on cell groups and a validated CSC/dense matrix.

Important internal types:

```rust
struct MetadataTable {
    row_names: Vec<String>,
    columns: Vec<MetadataColumn>,
}

struct CscMatrix {
    rows: usize,
    cols: usize,
    row_indices: Vec<u32>,
    column_offsets: Vec<u64>,
    values: Vec<f64>,
    feature_names: Vec<String>,
    cell_names: Vec<String>,
}
```

## Seurat interpretation rules

### Top-level object

Validate that the root is an S4 object inheriting from `Seurat`. Read only:

- `meta.data`
- `assays`
- `active.assay`
- names of `reductions` for inspection

Cell metadata row names are authoritative identifiers. Expression cells are
matched by name and reordered; positional equality alone is insufficient.

### Seurat v4

Read `assays[[assay]]@data` or the explicitly requested slot. Matrix `Dimnames`
provide features and cells.

### Seurat v5

Read `assays[[assay]]@layers[[layer]]`. A layer may contain only a subset of the
assay's cells and features. Reconstruct its names from the assay's `cells` and
`features` `LogMap` membership columns and honor the matrix orientation. Do not
assume every layer shares the full assay dimensions.

If multiple layers match `data` (for example `data.sample1` and `data.sample2`),
the first release reports ambiguity and asks for an exact layer. Joining layers
is a separate, tested feature.

## Sparse averaging

For log1p-normalized CSC input, allocate the unavoidable dense output buffer of
`features * groups` doubles and iterate columns once:

```text
for each cell column c:
    g = group_for_cell[c]
    group_size[g] += 1
    for each stored (gene, value) in c:
        sums[gene, g] += expm1(value)

for each gene and group:
    result[gene, g] = log1p(sums[gene, g] / group_size[g])
```

Implicit sparse zeroes contribute zero to the sum but their cells remain in the
divisor. Linear input skips `expm1`. Accumulation uses `f64`.

## Large-file strategy

The supplied real fixture is approximately 1.9 GB, so full-object
materialization is not acceptable as the default architecture.

1. Inspect the RDS stream lazily to discover object paths and vector metadata.
2. Materialize metadata and small structural slots only.
3. Materialize or stream only the selected matrix's `i`, `p`, and `x` vectors.
4. Never materialize reductions, graphs, images, command logs, or unused assays.
5. Apply parser limits for nesting, vector sizes, and requested output size.

If `rds2rust` cannot selectively extract vectors nested in the Seurat S4 graph,
extend or patch that layer before implementing the full CLI around an in-memory
parse.

## Milestones and acceptance criteria

### M0: compatibility probe

- Parse headers and lazily traverse the real v5 fixture.
- Locate `meta.data`, the active assay, layers, and matrix vector paths.
- Extract dimensions and names without loading the expression values.
- Record unsupported R object types rather than panicking.

Go/no-go: proceed only when paths and dimensions agree with R.

### M1: metadata

- Implement `inspect`, `metadata --columns`, and `head`.
- Match R's metadata values and row names on fixtures.
- Render missing values consistently as `NA`.

### M2: expression layer

- Validate `dgCMatrix` invariants (`p`, `i`, `x`, `Dim`, `Dimnames`).
- Reconstruct v4 and v5 layer names and orientation.
- Match selected entries and dimensions reported by R.

### M3: reference matrix

- Align metadata to expression cells by name.
- Define and test factor ordering, `NA` annotations, and empty levels.
- Match a pinned clustifyr version on small golden fixtures within a documented
  floating-point tolerance.
- Produce deterministic, atomic TSV output.

### M4: hardening and releases

- Corrupt/truncated RDS tests and allocation-limit tests.
- Benchmarks for elapsed time and peak resident memory.
- Linux x86_64/ARM64 and macOS x86_64/ARM64 release binaries.
- Checksums and smoke tests that run without R installed.

## Test matrix

- v4 and v5 Seurat objects
- factor and character annotations, including unused levels and `NA`
- reordered metadata rows
- partial v5 cell/feature layers
- multiple assays and non-RNA active assay
- single and split layers
- dense and `dgCMatrix` data
- duplicate/missing cell and feature names
- gzip, bzip2, xz, and uncompressed RDS envelopes as supported
- unsupported external-backed matrix with a clear diagnostic

Small fixtures and expected TSV files should be versioned. Large real-world RDS
files should remain outside Git or use an explicit large-data mechanism.

## Dependency policy

Start with `clap`, pinned `rds2rust`, `anyhow`/`thiserror`, and `csv`. Add `sprs`
only if it simplifies validated CSC access; grouped averaging itself does not
require it. Pin all dependencies in `Cargo.lock` and audit release licenses.

## Current proof target

Use `testdata/so.rds` to answer the central question: can the current
`rds2rust` APIs selectively find and extract Seurat v5 metadata and an ordinary
normalized expression layer at this file size? The answer to that question
determines whether the next change is the CLI implementation or a focused RDS
reader extension.
