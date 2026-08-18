# Refmat: Standalone Seurat Reference Matrix CLI

## Goal

Build a small, fast, self-contained CLI named `refmat` that can directly read R `.rds` files containing Seurat objects **without requiring R, Seurat, or clustifyr to be installed**.

The tool should initially provide two core capabilities:

1. Inspect a Seurat object and print metadata.
2. Generate a clustifyr-compatible reference expression matrix from a Seurat metadata annotation column.

Rust is the preferred implementation language so releases can be distributed as single binaries for Linux and macOS.

## Desired CLI

### Inspect metadata

```
refmat head so.rds
```

Equivalent conceptually to:

```
head(so@meta.data)
```

Support:

```
refmat head so.rds -n 20
```

Additional useful commands:

```
refmat inspect so.rds
refmat metadata so.rds --columns
```

`inspect` should summarize the object, for example:

```
Seurat object
Cells:     84,291
Features:  36,601

Active assay: RNA

Assays:
  RNA
    counts  36,601 × 84,291
    data    36,601 × 84,291

Reductions:
  pca
  harmony
  umap

Metadata:
  orig.ident
  nCount_RNA
  nFeature_RNA
  sample
  seurat_clusters
  cell_types
```

## Reference Matrix Generation

Primary command:

```
refmat build so.rds --column cell_types
```

Potential shorthand:

```
refmat so.rds --column cell_types
```

The command should generate something equivalent to the reference matrix produced by `clustifyr` from a Seurat object.

For each unique value in:

```
so@meta.data$cell_types
```

calculate the average expression profile for all cells belonging to that group.

For normalized Seurat expression data, clustifyr essentially performs:

```
mat_data <- expm1(mat[, cells])
res <- Matrix::rowMeans(mat_data)
res <- log1p(res)
```

Therefore the Rust implementation should reproduce:

```
normalized expression matrix
        ↓
group cells using metadata column
        ↓
for each group
        ↓
expm1(expression)
        ↓
row mean
        ↓
log1p(mean)
        ↓
genes × cell types reference matrix
```

Regression tests should compare the Rust output against the equivalent clustifyr output.

## RDS Parsing

Use the Rust crate:

```
rds2rust
```

as the initial RDS deserialization implementation.

The major reason this project is now practical is that `rds2rust` supports many R structures required by Seurat, including:

- S4 objects
- attributes
- lists
- environments
- matrices
- factors
- reference tracking
- ALTREP

Do **not** embed an R runtime unless pure-Rust parsing proves impossible.

The desired architecture is:

```
RDS
 ↓
rds2rust
 ↓
generic R object representation
 ↓
Seurat-specific parser
 ↓
internal Rust Seurat representation
```

## Seurat Parsing

The project does not need to reimplement Seurat.

It only needs to understand the subset of Seurat structures required for metadata and expression extraction.

Initially target:

```
Seurat v4
Seurat v5
```

Important structures include:

```
Seurat
├── meta.data
├── assays
│   └── RNA
│       ├── Assay (v4)
│       └── Assay5 (v5)
│           └── layers
│               ├── counts
│               └── data
├── active.assay
├── reductions
├── cell names
└── feature names
```

Create a clean internal abstraction rather than exposing raw R structures throughout the program.

Conceptually:

```
struct SeuratObject {
    metadata: Metadata,
    assays: HashMap<String, Assay>,
    active_assay: String,
}

struct Assay {
    features: Vec<String>,
    cells: Vec<String>,
    layers: HashMap<String, SparseMatrix>,
}
```

The exact structures should be designed based on what `rds2rust` exposes.

## Sparse Matrix Support

Seurat expression matrices are commonly `dgCMatrix` objects.

Important slots include:

```
i
p
x
Dim
Dimnames
```

These correspond naturally to CSC sparse matrices.

Use a Rust sparse matrix library such as:

```
sprs
```

and avoid converting the complete expression matrix to dense representation.

Reference averaging should operate directly against sparse data wherever possible.

This is important because real Seurat objects may contain tens or hundreds of thousands of cells.

## Assay and Layer Selection

Do not hard-code `RNA/data`.

Defaults should be:

```
assay = active assay
layer = data
```

Allow explicit overrides:

```
refmat build so.rds \
    --column cell_types \
    --assay RNA \
    --layer data
```

The program should produce informative errors when the requested assay, layer, or metadata column does not exist.

## Output

Start with a simple interoperable TSV representation:

```
refmat build so.rds --column cell_types
```

producing:

```
so.refmat.tsv
```

Example:

```
gene    CD4 T    CD8 T    NK    B cell    Monocyte
CD3D    5.12     5.34     0.12  0.03      0.04
CD3E    4.92     5.10     0.08  0.02      0.03
NKG7    1.21     4.81     6.22  0.10      1.94
```

Optional future formats could include:

```
Parquet
Arrow
RDS
```

TSV should be sufficient for the first implementation and makes validation against R straightforward.

Allow:

```
refmat build so.rds \
    --column cell_types \
    --output PBMC.refmat.tsv
```

## Suggested Rust Dependencies

Investigate:

```
clap        CLI
rds2rust    RDS deserialization
sprs        sparse matrices
ndarray     numerical operations if necessary
csv         TSV/CSV writing
anyhow      error handling
```

Avoid unnecessary dependencies.

## Development Milestones

### Milestone 1 — RDS/Seurat Explorer

Take real Seurat v4/v5 `.rds` files and determine whether `rds2rust` can successfully deserialize them.

Implement:

```
refmat inspect so.rds
```

Verify access to:

- Seurat S4 object
- `meta.data`
- assays
- `Assay` / `Assay5`
- layers
- `dgCMatrix`
- dimensions
- feature names
- cell names

This is the primary technical proof-of-concept.

### Milestone 2 — Metadata

Implement:

```
refmat head so.rds
refmat head so.rds -n 20
refmat metadata so.rds --columns
```

Correctly handle common R dataframe column types:

- character
- factor
- integer
- numeric
- logical
- missing values

### Milestone 3 — Expression Matrix Extraction

Implement internal extraction of:

```
--assay RNA
--layer data
```

Correctly reconstruct `dgCMatrix` as a Rust CSC sparse matrix.

Verify:

```
dimensions
feature names
cell names
values
```

against R.

### Milestone 4 — Reference Matrix

Implement:

```
refmat build so.rds --column cell_types
```

Perform grouped expression averaging equivalent to clustifyr.

Validate using multiple real Seurat objects.

For example:

```
rust <- read.delim("so.refmat.tsv", row.names = 1)

clustifyr_ref <- clustifyr::seurat_ref(
    so,
    cluster_col = "cell_types"
)

all.equal(
    as.matrix(rust),
    as.matrix(clustifyr_ref),
    tolerance = 1e-10
)
```

The exact clustifyr API and calculation should be verified against the current clustifyr source during implementation.

## Testing Strategy

Create small Seurat test objects in R covering:

```
Seurat v4
Seurat v5
multiple assays
Assay5 layers
factor annotations
character annotations
NA annotations
empty groups
different active assays
sparse matrices
large sparse matrices
```

Store expected metadata and reference matrices alongside the test fixtures.

Rust integration tests should compare against those known outputs.

Real-world Seurat objects should also be tested before considering the parser stable.

## Distribution

The final goal is a standalone executable:

```
refmat
```

with GitHub Releases containing at least:

```
Linux x86_64
Linux ARM64
macOS ARM64
macOS x86_64
```

Users should be able to:

```
chmod +x refmat
./refmat inspect so.rds
```

without installing:

```
R
Seurat
clustifyr
Python
```

## Important Design Principle

Keep these layers separated:

```
RDS parsing
     ↓
Seurat interpretation
     ↓
expression/metadata abstraction
     ↓
reference matrix algorithm
     ↓
CLI/output
```

Do not intermingle `rds2rust` internals with the reference-matrix algorithm.

Ideally the reference matrix code operates on generic structures such as:

```
metadata groups
feature names
cell names
CSC expression matrix
```

This will make testing much easier and potentially allow support for other single-cell formats later.

## Primary Technical Risk

The biggest risk is not reference-matrix generation.

It is correctly interpreting complex Seurat S4 structures across Seurat versions.

Therefore **do not begin by implementing the entire CLI**.

The first prototype should answer one question:

> Can `rds2rust` take a real Seurat v5 `.rds` file and reliably extract `meta.data` and reconstruct the normalized expression `dgCMatrix` from an `Assay5` layer?

If yes, proceed with the full Rust implementation.

If not, investigate extending/fixing `rds2rust` before considering embedding R or switching implementation strategies.

## Initial Target

The first working version only needs to make these commands reliable:

```
refmat inspect so.rds

refmat head so.rds

refmat build so.rds --column cell_types
```

Everything else should be considered secondary until these work correctly against real Seurat objects.