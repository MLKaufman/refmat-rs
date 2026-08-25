# refmat


<p align="center">
<img src="logo.png" alt="refmat logo" width="50%"/>
</p>

`refmat` is a standalone command-line tool for inspecting single-cell objects
and building genes-by-cell-type reference expression matrices. 

<center><b> Reads Seurat, SingleCellExperiment, and AnnData files directly, without requiring R or Python
at runtime. </center></b>

## Why

`refmat` was created to provide a fast, efficient way to build cell-type reference expression matrices from single-cell data without the overhead of running R or Python scripts. This makes it ideal for use in automated pipelines and for users who want to avoid the complexity of managing and maintaining R or Python environments.

**Want to explore existing single-cell dataset objects, quickly summarize their contents, and generate reference matrices for downstream cell type annotation analysis?**

This  tool to makes it easy to do all of that from the command line, with minimal dependencies and maximum speed.

To help support this workflow during single-cell analysis:  
**Single Cell Experiment Object -> Reference Matrix -> ClustifyR -> Cell Type Annotation**

Please check out the [ClustifyR](https://github.com/rnabioco/clustifyr/) or its python port [pyclustifyr](https://github.com/MLKaufman/pyclustifyr) repositories for more information on using a downstream reference matrix for cell type annotation.

## Highlights

- Automatically detects Seurat and SingleCellExperiment RDS files and AnnData
  H5AD files from their contents.
- Previews cell metadata in readable, terminal-sized tables.
- Counts cells by a categorical or character metadata column.
- Aggregates dense and sparse expression matrices with bounded memory use.
- Produces a clustifyr-compatible genes-by-group TSV reference matrix.

## Supported inputs

| Input | Cell metadata | Expression matrices |
| --- | --- | --- |
| Seurat v3/v4 RDS | `meta.data` | Legacy `Assay` sparse or dense `counts`, `data`, and `scale.data` slots |
| Seurat v5 RDS | `meta.data` | In-memory `Assay5` `dgCMatrix` layers |
| SingleCellExperiment RDS | `colData` | Dense integer/double matrices and `dgCMatrix` assays |
| AnnData H5AD | `obs` | Dense, CSR, or CSC `X` and named layers, with gzip or LZF compression |

*File type detection is content-based; the filename extension does not determine
whether an input is treated as H5AD or RDS.*

## Installation

### Source via Cargo

Install the latest release directly from GitHub:

```bash
cargo install \
  --git https://github.com/MLKaufman/refmat-rs.git \
  --tag v1.1.1 \
  --locked
```

*Note: Building requires Rust 1.85 or newer, CMake, and a C compiler. HDF5, gzip/LZF,
and the RDS xz decoder are built into the executable, so running the installed
binary does not require R, Python, a system HDF5 installation, or liblzma.*

### Prebuilt binaries

Download the latest release for your platform from the [GitHub releases page](https://github.com/MLKaufman/refmat-rs/releases).

Prebuilt binaries are available for Linux, and ARM-based macOS. After downloading, add the binary to your PATH.

## Quick start

```bash
# Summarize an object.
refmat inspect so.rds

# Preview six rows of cell metadata.
refmat head so.rds

# Count cells by annotation.
refmat col so.rds cell_type

# Build sample.refmat.tsv from that annotation.
refmat build so.rds cell_type
```

Run `refmat --help` or `refmat <COMMAND> --help` for the complete CLI help.

## Commands

### Inspect an object

```bash
refmat inspect <FILE>
refmat inspect <FILE> --depth 6
```

`inspect` summarizes the serialized object without loading its expression
matrix into memory. For small diagnostic RDS files, `--full` materializes all
vectors. Seurat inputs report their stored object version, active assay, and
whether that assay uses the v3/v4 `Assay` layout or the v5 `Assay5` layout.

### Preview cell metadata

```bash
refmat head <FILE>
refmat head <FILE> --rows 20
```

The default is six cells. Small metadata frames print as one aligned table;
wide frames are split into labeled table blocks that repeat the cell identifier
and remain readable in a terminal. For Seurat inputs, a detection line reporting
the stored version and active-assay layout is written to stderr.

### Count cells by metadata value

```bash
refmat col <FILE> <COLUMN>
refmat col <FILE> --column <COLUMN>
```

For example:

```bash
refmat col so.rds cell_type
```

The column may be positional or supplied with `-c`/`--column`; for example,
`refmat col so.rds -c cell_type` is equivalent.

```text
+-----------+-------+
| cell_type | cells |
+-----------+-------+
| B cell    | 2     |
| T cell    | 2     |
+-----------+-------+
```

The grouping column must be a factor or character vector in RDS, or a
categorical or string column in H5AD. Unused factor/category levels are shown
with a count of zero; missing annotations are excluded.

### Build a reference matrix

```bash
refmat build <FILE> <COLUMN>
refmat build <FILE> --column <COLUMN>
```

As with `col`, the column may be positional or supplied with `-c`/`--column`.

Defaults depend on the input format:

| Input | Default expression matrix | Overrides |
| --- | --- | --- |
| Seurat | Active assay, `data` layer | `--assay`, `--layer` |
| SingleCellExperiment | `logcounts`, then `counts`, then the first assay | `--assay` |
| AnnData | `X` | `--layer` |

Select inputs and output explicitly when needed:

```bash
# Seurat
refmat build so.rds \
  --column celltypes \
  --assay RNA \
  --layer data \
  --output custom.refmat.tsv

# SingleCellExperiment
refmat build sce.rds \
  --column cell_type \
  --assay logcounts

# AnnData
refmat build sample.h5ad \
  --column cell_type \
  --layer counts \
  --scale linear
```

Without `--output`, the result is written beside the input as
`<input-stem>.refmat.tsv`. Rows are features, columns are annotation groups,
and the first column is named `gene`.

For Seurat v3/v4 `Assay` objects, `--layer` selects the corresponding direct
slot and must be `counts`, `data`, or `scale.data`. For Seurat v5 `Assay5`
objects, it selects a named layer.

## Expression scale and aggregation

For each feature and annotation group, `refmat` calculates:

```text
log1p(mean(linear expression))
```

Log1p-normalized input is converted with `expm1` before averaging. Linear input,
such as raw counts, is averaged directly before `log1p` is applied.

The default `--scale auto` treats a matrix named `counts` as linear and all
other assay/layer names as log1p-normalized. AnnData does not record whether
`X` contains counts or normalized values, so specify the scale when necessary:

```bash
refmat build raw-counts.h5ad --column cell_type --scale linear
refmat build normalized.h5ad --column cell_type --scale log1p
```

Character annotations preserve first-seen order. R factors preserve factor
level order, and AnnData categoricals preserve category order. Cells with
missing annotations are excluded from the reference matrix.

## Validation

The test suite includes structurally faithful Seurat v3/v4 fixtures and genuine
SingleCellExperiment and AnnData fixtures. Together they exercise version and
assay-layout detection, metadata decoding, dense matrices, sparse matrices,
categorical group ordering, and automatic scale handling.

The Seurat implementation was additionally validated on a 61,844-cell,
32,285-feature Seurat v5 object with 135,265,589 nonzero entries. The generated
32,285 × 18 reference matrix matched the equivalent R/Matrix calculation:

```text
maximum absolute error:  8.881784e-16
mean absolute error:     2.137444e-18
all.equal tolerance 1e-12: TRUE
```

The large validation object is not included in the repository. The comparison
script is available at `scripts/validate_reference.R` for compatible Seurat
objects.

## Current limitations

- Seurat v5 `Assay5` expression layers must be `dgCMatrix`. Legacy v3/v4
  `Assay` slots may be `dgCMatrix` or dense integer/double matrices.
- Seurat layers must contain all cells in metadata order. Partial and split
  layers are not aligned or joined automatically.
- BPCells, DelayedArray/HDF5Array within RDS, and external-pointer matrix
  backends are not supported. Export these objects to an in-memory assay or
  H5AD first.
- AnnData `.raw/X` is not selectable; use `X` or a named layer.
- H5AD grouping columns must be categorical or string values.
- Each invocation decompresses a compressed RDS into a temporary backing file;
  repeated commands do not currently share a persistent cache.

## License

`refmat` is licensed under the MIT License.
