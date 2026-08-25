use std::collections::HashMap;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, anyhow, bail, ensure};
use clap::{Parser, Subcommand, ValueEnum};
use rds2rust::{
    Attributes, ChunkedRdsSource, Logical, ParseConfig, RObject, RdsInput, VectorData,
    read_lazy_character_range, read_lazy_integer_range, read_lazy_logical_range,
    read_lazy_real_range, read_rds_from_path_chunked, read_rds_with_input,
};

mod h5ad;

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum InputScale {
    /// Infer from the assay/layer name; otherwise assume log1p-normalized values.
    Auto,
    /// Values are log1p-normalized and are averaged on the linear scale.
    Log1p,
    /// Values are already linear (for example raw counts).
    Linear,
}

#[derive(Debug, Parser)]
#[command(version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Summarize the serialized object without loading large vectors.
    Inspect {
        file: PathBuf,
        /// Print the generic R object tree to this depth.
        #[arg(long, default_value_t = 4)]
        depth: usize,
        /// Materialize all vectors. Intended only for small diagnostic files.
        #[arg(long)]
        full: bool,
    },
    /// Print the first rows of cell metadata.
    Head {
        file: PathBuf,
        #[arg(short = 'n', long, default_value_t = 6)]
        rows: usize,
    },
    /// Count cells by the values in a metadata column.
    Col {
        file: PathBuf,
        /// Metadata column (may also be supplied with -c/--column).
        #[arg(value_name = "COLUMN", conflicts_with = "column_option")]
        column: Option<String>,
        /// Metadata column (alternative to the positional COLUMN).
        #[arg(short = 'c', long = "column", value_name = "COLUMN")]
        column_option: Option<String>,
    },
    /// Build a genes-by-group reference matrix.
    Build {
        file: PathBuf,
        /// Metadata column (may also be supplied with -c/--column).
        #[arg(value_name = "COLUMN", conflicts_with = "column_option")]
        column: Option<String>,
        /// Metadata column (alternative to the positional COLUMN).
        #[arg(short = 'c', long = "column", value_name = "COLUMN")]
        column_option: Option<String>,
        /// Seurat or SingleCellExperiment assay name.
        #[arg(long)]
        assay: Option<String>,
        /// Seurat or H5AD layer; also selects a legacy Assay matrix slot.
        #[arg(long)]
        layer: Option<String>,
        /// Interpretation of matrix values before group averaging.
        #[arg(long, value_enum, default_value_t = InputScale::Auto)]
        scale: InputScale,
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Inspect { file, depth, full } => inspect(&file, depth, full),
        Command::Head { file, rows } => head(&file, rows),
        Command::Col {
            file,
            column,
            column_option,
        } => col(&file, &required_column(column, column_option)?),
        Command::Build {
            file,
            column,
            column_option,
            assay,
            layer,
            scale,
            output,
        } => build(
            &file,
            &required_column(column, column_option)?,
            assay.as_deref(),
            layer.as_deref(),
            scale,
            output.as_deref(),
        ),
    }
}

fn required_column(positional: Option<String>, option: Option<String>) -> Result<String> {
    positional.or(option).ok_or_else(|| {
        anyhow!("a metadata column is required; provide COLUMN or -c/--column COLUMN")
    })
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum InputFormat {
    Rds,
    H5ad,
}

fn detect_format(file: &Path) -> Result<InputFormat> {
    let mut handle =
        File::open(file).with_context(|| format!("failed to open {}", file.display()))?;
    let mut magic = [0u8; 8];
    let read = handle.read(&mut magic)?;
    if read == magic.len() && magic == [0x89, b'H', b'D', b'F', b'\r', b'\n', 0x1a, b'\n'] {
        Ok(InputFormat::H5ad)
    } else {
        Ok(InputFormat::Rds)
    }
}

struct ParsedRds {
    // The source owns the temporary decompressed file that lazy vector spans
    // refer to, so it must outlive every operation on `object`.
    source: ChunkedRdsSource,
    object: RObject,
    lazy_vectors: usize,
}

fn parse_lazy(file: &Path) -> Result<ParsedRds> {
    // Seurat cell metadata commonly has dozens of columns. The crate's generic
    // lazy default skips list contents longer than ten, which preserves the
    // row count but loses all dataframe columns. The trusted-large preset keeps
    // lists up to 100 elements structural while leaving their long vectors lazy.
    let source = ChunkedRdsSource::from_path(file)
        .with_context(|| format!("failed to open {}", file.display()))?;
    let parsed = read_rds_with_input(&source, ParseConfig::for_trusted_large_file())
        .with_context(|| format!("failed to parse {}", file.display()))?;
    Ok(ParsedRds {
        source,
        object: parsed.object,
        lazy_vectors: parsed.warnings.len(),
    })
}

fn inspect(file: &Path, depth: usize, full: bool) -> Result<()> {
    if detect_format(file)? == InputFormat::H5ad {
        return h5ad::inspect(file, depth);
    }
    if full {
        let object = read_rds_from_path_chunked(file)
            .with_context(|| format!("failed to parse {}", file.display()))?
            .object;
        if is_seurat(&object) {
            println!("{}", seurat_format_summary(&object)?);
        }
        print_object("root", &object, 0, depth);
    } else {
        let parsed = parse_lazy(file)?;
        eprintln!("note: {} large vectors kept lazy", parsed.lazy_vectors);
        if is_seurat(&parsed.object) {
            println!("{}", seurat_format_summary(&parsed.object)?);
        }
        print_object("root", &parsed.object, 0, depth);
    }
    Ok(())
}

fn head(file: &Path, rows: usize) -> Result<()> {
    if detect_format(file)? == InputFormat::H5ad {
        return h5ad::head(file, rows);
    }
    let parsed = parse_lazy(file)?;
    if is_sce(&parsed.object) {
        return sce_head(&parsed, rows);
    }
    eprintln!("{}", seurat_format_summary(&parsed.object)?);
    let frame = seurat_metadata_frame(&parsed.object)?;
    let take = rows.min(frame.row_names.len());
    let row_names = seurat_cell_names(&parsed.object, &parsed.source, take)
        .unwrap_or_else(|_| frame.row_names.iter().take(take).cloned().collect());

    let rendered = frame
        .columns
        .iter()
        .map(|(name, column)| {
            render_column_range(column, &parsed.source, take)
                .with_context(|| format!("failed to read metadata column '{name}'"))
        })
        .collect::<Result<Vec<_>>>()?;

    let row_names = row_names
        .iter()
        .map(|name| name.as_deref().unwrap_or("NA").to_owned())
        .collect::<Vec<_>>();
    print_metadata_table(
        &row_names,
        frame.columns.keys().map(|name| name.as_ref()),
        &rendered,
    );
    Ok(())
}

fn col(file: &Path, column: &str) -> Result<()> {
    if detect_format(file)? == InputFormat::H5ad {
        return h5ad::col(file, column);
    }
    let parsed = parse_lazy(file)?;
    let annotation = if is_sce(&parsed.object) {
        sce_columns(&parsed.object)?
            .into_iter()
            .find_map(|(name, values)| (name == column).then_some(values))
            .ok_or_else(|| anyhow!("colData column '{column}' does not exist"))?
    } else {
        seurat_metadata_frame(&parsed.object)?
            .columns
            .get(column)
            .ok_or_else(|| anyhow!("metadata column '{column}' does not exist"))?
    };
    let (group_names, cell_groups) = read_groups(annotation, &parsed.source)?;
    ensure!(
        !group_names.is_empty(),
        "metadata column '{column}' has no non-missing groups"
    );
    print_group_counts(column, &group_names, &cell_groups);
    Ok(())
}

fn build(
    file: &Path,
    column: &str,
    assay: Option<&str>,
    layer: Option<&str>,
    scale: InputScale,
    output: Option<&Path>,
) -> Result<()> {
    if detect_format(file)? == InputFormat::H5ad {
        ensure!(assay.is_none(), "--assay is only valid for RDS inputs");
        return h5ad::build(file, column, layer, scale, output);
    }
    let parsed = parse_lazy(file)?;
    if is_sce(&parsed.object) {
        ensure!(
            layer.is_none(),
            "--layer is only valid for Seurat and H5AD inputs"
        );
        return sce_build(&parsed, file, column, assay, scale, output);
    }
    let frame = seurat_metadata_frame(&parsed.object)?;
    let annotation = frame
        .columns
        .get(column)
        .ok_or_else(|| anyhow!("metadata column '{column}' does not exist"))?;
    let (group_names, cell_groups) = read_groups(annotation, &parsed.source)?;
    ensure!(
        !group_names.is_empty(),
        "metadata column '{column}' has no non-missing groups"
    );

    let active_assay = active_assay_name(&parsed.object)?;
    let assay_name = assay.unwrap_or(&active_assay);
    let layer_name = layer.unwrap_or("data");
    let assay_s4 = seurat_assay(&parsed.object, assay_name)?;
    let layout = seurat_assay_layout(assay_s4)?;
    let (matrix_object, feature_names, cell_names) = match layout {
        SeuratAssayLayout::V5 => {
            let matrix = named_item(
                assay_s4
                    .slots
                    .get("layers")
                    .ok_or_else(|| anyhow!("Assay5 has no layers slot"))?,
                layer_name,
            )
            .with_context(|| {
                format!("layer '{layer_name}' does not exist in assay '{assay_name}'")
            })?;
            let cells = layer_member_names(assay_s4, "cells", layer_name, &parsed.source)?;
            let features = layer_member_names(assay_s4, "features", layer_name, &parsed.source)?;
            (matrix, features, cells)
        }
        SeuratAssayLayout::Legacy => {
            ensure!(
                matches!(layer_name, "counts" | "data" | "scale.data"),
                "legacy Assay layer must be one of counts, data, or scale.data"
            );
            let matrix = assay_s4
                .slots
                .get(layer_name)
                .ok_or_else(|| anyhow!("legacy assay '{assay_name}' has no '{layer_name}' slot"))?;
            let features = r_matrix_names(matrix, &parsed.source, 0)
                .with_context(|| format!("failed to read feature names from '{layer_name}'"))?;
            let cells = r_matrix_names(matrix, &parsed.source, 1)
                .with_context(|| format!("failed to read cell names from '{layer_name}'"))?;
            (matrix, features, cells)
        }
    };

    ensure_metadata_cell_order(frame, &cell_names, assay_name, layer_name)?;
    let sums = match sparse_matrix(matrix_object) {
        Ok(matrix) => {
            validate_matrix_names(&matrix, &feature_names, &cell_names, layer_name)?;
            aggregate_r_sparse(
                &matrix,
                &parsed.source,
                &cell_groups,
                group_names.len(),
                scale,
                layer_name,
            )?
        }
        Err(sparse_error) if layout == SeuratAssayLayout::Legacy => {
            let matrix = dense_matrix(matrix_object).with_context(|| {
                format!(
                    "legacy assay '{assay_name}' slot '{layer_name}' is neither dgCMatrix nor a dense numeric matrix ({sparse_error})"
                )
            })?;
            validate_matrix_names(&matrix, &feature_names, &cell_names, layer_name)?;
            aggregate_r_dense(
                &matrix,
                &parsed.source,
                &cell_groups,
                group_names.len(),
                scale,
                layer_name,
            )?
        }
        Err(error) => return Err(error),
    };

    let output_path = output
        .map(Path::to_path_buf)
        .unwrap_or_else(|| default_output(file));
    write_reference(&output_path, &feature_names, &group_names, &sums)?;
    eprintln!(
        "wrote {} features x {} groups to {}",
        feature_names.len(),
        group_names.len(),
        output_path.display()
    );
    Ok(())
}

struct SparseMatrixRef<'a> {
    rows: usize,
    cols: usize,
    i: &'a VectorData<i32>,
    p: &'a VectorData<i32>,
    x: &'a VectorData<f64>,
}

enum DenseValuesRef<'a> {
    Real(&'a VectorData<f64>),
    Integer(&'a VectorData<i32>),
}

struct DenseMatrixRef<'a> {
    rows: usize,
    cols: usize,
    values: DenseValuesRef<'a>,
}

trait MatrixDimensions {
    fn rows(&self) -> usize;
    fn cols(&self) -> usize;
}

impl MatrixDimensions for SparseMatrixRef<'_> {
    fn rows(&self) -> usize {
        self.rows
    }

    fn cols(&self) -> usize {
        self.cols
    }
}

impl MatrixDimensions for DenseMatrixRef<'_> {
    fn rows(&self) -> usize {
        self.rows
    }

    fn cols(&self) -> usize {
        self.cols
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SeuratAssayLayout {
    Legacy,
    V5,
}

impl SeuratAssayLayout {
    fn description(self) -> &'static str {
        match self {
            Self::Legacy => "Assay (v3/v4 layout)",
            Self::V5 => "Assay5 (v5 layout)",
        }
    }
}

fn is_seurat(object: &RObject) -> bool {
    as_s4(object, "RDS root")
        .is_ok_and(|s4| s4.class.iter().any(|class| class.as_ref() == "Seurat"))
}

fn is_sce(object: &RObject) -> bool {
    as_s4(object, "RDS root").is_ok_and(|s4| {
        s4.class
            .iter()
            .any(|class| class.as_ref() == "SingleCellExperiment")
    })
}

fn sce_coldata(object: &RObject) -> Result<&rds2rust::S4ObjectData> {
    let root = as_s4(object, "SingleCellExperiment")?;
    ensure!(
        is_sce(object),
        "RDS object is neither Seurat nor SingleCellExperiment"
    );
    as_s4(
        root.slots
            .get("colData")
            .ok_or_else(|| anyhow!("SingleCellExperiment has no colData slot"))?,
        "SingleCellExperiment colData",
    )
}

fn sce_columns(object: &RObject) -> Result<Vec<(String, &RObject)>> {
    let frame = sce_coldata(object)?;
    let list = frame
        .slots
        .get("listData")
        .ok_or_else(|| anyhow!("SingleCellExperiment colData has no listData slot"))?;
    let (value, attrs) = with_attributes(list)?;
    let values = match value {
        RObject::List(values) => values,
        other => bail!(
            "SingleCellExperiment colData listData decoded as {}",
            describe(other)
        ),
    };
    let names = loaded_characters(
        attrs
            .get("names")
            .ok_or_else(|| anyhow!("colData has no column names"))?,
        "colData names",
    )?;
    ensure!(
        names.len() == values.len(),
        "colData names and values have different lengths"
    );
    Ok(names
        .iter()
        .zip(values)
        .map(|(name, value)| (name.as_deref().unwrap_or("NA").to_owned(), value))
        .collect())
}

fn sce_row_names(
    object: &RObject,
    source: &dyn RdsInput,
    count: usize,
) -> Result<Vec<Option<Arc<str>>>> {
    let frame = sce_coldata(object)?;
    let names = frame
        .slots
        .get("rownames")
        .ok_or_else(|| anyhow!("SingleCellExperiment colData has no rownames"))?;
    read_character_range(names, source, 0, count)
}

fn sce_nrows(object: &RObject) -> Result<usize> {
    let frame = sce_coldata(object)?;
    let values = loaded_integers(
        frame
            .slots
            .get("nrows")
            .ok_or_else(|| anyhow!("colData has no nrows slot"))?,
        "colData nrows",
    )?;
    let value = *values
        .first()
        .ok_or_else(|| anyhow!("colData nrows is empty"))?;
    ensure!(value >= 0, "colData nrows is negative");
    Ok(value as usize)
}

fn sce_head(parsed: &ParsedRds, rows: usize) -> Result<()> {
    let columns = sce_columns(&parsed.object)?;
    let take = rows.min(sce_nrows(&parsed.object)?);
    let row_names = sce_row_names(&parsed.object, &parsed.source, take)
        .unwrap_or_else(|_| (1..=take).map(|i| Some(Arc::from(i.to_string()))).collect());
    let rendered = columns
        .iter()
        .map(|(name, column)| {
            render_column_range(column, &parsed.source, take)
                .with_context(|| format!("failed to read colData column '{name}'"))
        })
        .collect::<Result<Vec<_>>>()?;

    let row_names = row_names
        .iter()
        .map(|name| name.as_deref().unwrap_or("NA").to_owned())
        .collect::<Vec<_>>();
    print_metadata_table(
        &row_names,
        columns.iter().map(|(name, _)| name.as_str()),
        &rendered,
    );
    Ok(())
}

fn sce_assays(object: &RObject) -> Result<&RObject> {
    let root = as_s4(object, "SingleCellExperiment")?;
    let assays = as_s4(
        root.slots
            .get("assays")
            .ok_or_else(|| anyhow!("SingleCellExperiment has no assays slot"))?,
        "SingleCellExperiment assays",
    )?;
    let data = as_s4(
        assays
            .slots
            .get("data")
            .ok_or_else(|| anyhow!("SimpleAssays has no data slot"))?,
        "SingleCellExperiment assay list",
    )?;
    data.slots
        .get("listData")
        .ok_or_else(|| anyhow!("assay list has no listData slot"))
}

fn sce_assay_names(object: &RObject) -> Result<Vec<String>> {
    let (_, attrs) = with_attributes(sce_assays(object)?)?;
    Ok(loaded_characters(
        attrs
            .get("names")
            .ok_or_else(|| anyhow!("SingleCellExperiment assays have no names"))?,
        "assay names",
    )?
    .iter()
    .map(|x| x.as_deref().unwrap_or("NA").to_owned())
    .collect())
}

fn r_matrix_names(matrix: &RObject, source: &dyn RdsInput, axis: usize) -> Result<Vec<String>> {
    let dimnames = if let Ok(s4) = as_s4(matrix, "matrix") {
        match s4.slots.get("Dimnames") {
            Some(RObject::List(values)) => values.as_slice(),
            _ => bail!("R matrix has no Dimnames"),
        }
    } else {
        let (_, attrs) = with_attributes(matrix)?;
        match attrs.get("dimnames") {
            Some(RObject::List(values)) => values.as_slice(),
            _ => bail!("R matrix has no dimnames attribute"),
        }
    };
    let values = dimnames
        .get(axis)
        .ok_or_else(|| anyhow!("R matrix Dimnames is incomplete"))?;
    Ok(read_all_characters(values, source)?
        .into_iter()
        .map(|x| x.as_deref().unwrap_or("NA").to_owned())
        .collect())
}

fn validate_matrix_names(
    matrix: &impl MatrixDimensions,
    feature_names: &[String],
    cell_names: &[String],
    matrix_name: &str,
) -> Result<()> {
    ensure!(
        matrix.rows() == feature_names.len(),
        "matrix '{matrix_name}' has {} rows but {} feature names",
        matrix.rows(),
        feature_names.len()
    );
    ensure!(
        matrix.cols() == cell_names.len(),
        "matrix '{matrix_name}' has {} columns but {} cell names",
        matrix.cols(),
        cell_names.len()
    );
    Ok(())
}

fn ensure_metadata_cell_order(
    frame: &rds2rust::DataFrameData,
    cell_names: &[String],
    assay_name: &str,
    layer_name: &str,
) -> Result<()> {
    ensure!(
        frame.row_names.len() == cell_names.len(),
        "metadata has {} cells but assay '{assay_name}' matrix '{layer_name}' has {}; partial-layer cell alignment is not implemented yet",
        frame.row_names.len(),
        cell_names.len()
    );
    let compact_row_names = frame.row_names.iter().enumerate().all(|(index, name)| {
        name.as_deref()
            .is_some_and(|name| name.parse::<usize>() == Ok(index + 1))
    });
    if compact_row_names {
        return Ok(());
    }
    for (index, (metadata_name, matrix_name)) in frame.row_names.iter().zip(cell_names).enumerate()
    {
        let metadata_name = metadata_name
            .as_deref()
            .ok_or_else(|| anyhow!("metadata cell name {} is missing", index + 1))?;
        ensure!(
            metadata_name == matrix_name,
            "cell order mismatch at position {}: metadata has '{metadata_name}' but assay '{assay_name}' matrix '{layer_name}' has '{matrix_name}'",
            index + 1
        );
    }
    Ok(())
}

fn sce_build(
    parsed: &ParsedRds,
    file: &Path,
    column: &str,
    assay: Option<&str>,
    scale: InputScale,
    output: Option<&Path>,
) -> Result<()> {
    let columns = sce_columns(&parsed.object)?;
    let annotation = columns
        .iter()
        .find(|(name, _)| name == column)
        .map(|(_, value)| *value)
        .ok_or_else(|| anyhow!("colData column '{column}' does not exist"))?;
    let (group_names, cell_groups) = read_groups(annotation, &parsed.source)?;
    ensure!(
        !group_names.is_empty(),
        "colData column '{column}' has no non-missing groups"
    );

    let assay_names = sce_assay_names(&parsed.object)?;
    let default_assay = if assay_names.iter().any(|name| name == "logcounts") {
        "logcounts"
    } else if assay_names.iter().any(|name| name == "counts") {
        "counts"
    } else {
        assay_names
            .first()
            .map(String::as_str)
            .ok_or_else(|| anyhow!("SingleCellExperiment has no assays"))?
    };
    let assay_name = assay.unwrap_or(default_assay);
    let matrix_object = named_item(sce_assays(&parsed.object)?, assay_name)
        .with_context(|| format!("assay '{assay_name}' does not exist"))?;
    let feature_names = r_matrix_names(matrix_object, &parsed.source, 0)?;
    let sums = if let Ok(matrix) = sparse_matrix(matrix_object) {
        ensure!(
            matrix.cols == cell_groups.len(),
            "assay has {} cells but colData has {}",
            matrix.cols,
            cell_groups.len()
        );
        ensure!(
            feature_names.len() == matrix.rows,
            "assay feature names do not match matrix rows"
        );
        aggregate_r_sparse(
            &matrix,
            &parsed.source,
            &cell_groups,
            group_names.len(),
            scale,
            assay_name,
        )?
    } else {
        let matrix = dense_matrix(matrix_object)?;
        ensure!(
            matrix.cols == cell_groups.len(),
            "assay has {} cells but colData has {}",
            matrix.cols,
            cell_groups.len()
        );
        ensure!(
            feature_names.len() == matrix.rows,
            "assay feature names do not match matrix rows"
        );
        aggregate_r_dense(
            &matrix,
            &parsed.source,
            &cell_groups,
            group_names.len(),
            scale,
            assay_name,
        )?
    };

    let output_path = output
        .map(Path::to_path_buf)
        .unwrap_or_else(|| default_output(file));
    write_reference(&output_path, &feature_names, &group_names, &sums)?;
    eprintln!(
        "wrote {} features x {} groups to {}",
        feature_names.len(),
        group_names.len(),
        output_path.display()
    );
    Ok(())
}

fn dense_matrix(object: &RObject) -> Result<DenseMatrixRef<'_>> {
    let (value, attrs) = with_attributes(object)?;
    let dims = loaded_integers(
        attrs
            .get("dim")
            .ok_or_else(|| anyhow!("dense matrix has no dim attribute"))?,
        "dense matrix dimensions",
    )?;
    ensure!(
        dims.len() == 2 && dims[0] >= 0 && dims[1] >= 0,
        "invalid dense matrix dimensions"
    );
    let values = match value {
        RObject::Real(values) => DenseValuesRef::Real(values),
        RObject::Integer(values) => DenseValuesRef::Integer(values),
        other => bail!("dense assay decoded as {}, not numeric", describe(other)),
    };
    let rows = dims[0] as usize;
    let cols = dims[1] as usize;
    let len = match values {
        DenseValuesRef::Real(values) => values.len(),
        DenseValuesRef::Integer(values) => values.len(),
    };
    ensure!(
        len == rows * cols,
        "dense matrix dimensions do not match its values"
    );
    Ok(DenseMatrixRef { rows, cols, values })
}

fn aggregate_r_dense(
    matrix: &DenseMatrixRef<'_>,
    source: &dyn RdsInput,
    cell_groups: &[Option<usize>],
    group_count: usize,
    scale: InputScale,
    matrix_name: &str,
) -> Result<Vec<f64>> {
    let mut group_sizes = vec![0usize; group_count];
    for group in cell_groups.iter().flatten() {
        group_sizes[*group] += 1;
    }
    ensure!(
        group_sizes.iter().all(|size| *size > 0),
        "annotation contains an empty factor level"
    );
    let mut sums = vec![0.0; matrix.rows * group_count];
    for (cell, group) in cell_groups.iter().copied().enumerate() {
        let Some(group) = group else { continue };
        match matrix.values {
            DenseValuesRef::Real(values) => {
                for (feature, value) in
                    read_real_range(values, source, cell * matrix.rows, matrix.rows)?
                        .into_iter()
                        .enumerate()
                {
                    sums[feature * group_count + group] += to_linear(value, scale, matrix_name);
                }
            }
            DenseValuesRef::Integer(values) => {
                for (feature, value) in
                    read_integer_range(values, source, cell * matrix.rows, matrix.rows)?
                        .into_iter()
                        .enumerate()
                {
                    ensure!(value != i32::MIN, "dense integer assay contains NA");
                    sums[feature * group_count + group] +=
                        to_linear(value as f64, scale, matrix_name);
                }
            }
        }
    }
    finish_group_means(&mut sums, &group_sizes);
    Ok(sums)
}

fn aggregate_r_sparse(
    matrix: &SparseMatrixRef<'_>,
    source: &dyn RdsInput,
    cell_groups: &[Option<usize>],
    group_count: usize,
    scale: InputScale,
    matrix_name: &str,
) -> Result<Vec<f64>> {
    let mut group_sizes = vec![0usize; group_count];
    for group in cell_groups.iter().flatten() {
        group_sizes[*group] += 1;
    }
    ensure!(
        group_sizes.iter().all(|size| *size > 0),
        "annotation contains an empty factor level"
    );
    let p = read_all_integers(matrix.p, source)?;
    ensure!(p.len() == matrix.cols + 1, "invalid dgCMatrix p length");
    ensure!(
        p.first() == Some(&0) && p.windows(2).all(|x| x[0] <= x[1]),
        "invalid dgCMatrix column pointers"
    );
    let nnz = *p.last().unwrap_or(&0);
    ensure!(
        nnz >= 0 && nnz as usize == matrix.i.len() && nnz as usize == matrix.x.len(),
        "invalid dgCMatrix nonzero lengths"
    );

    let mut sums = vec![0.0f64; matrix.rows * group_count];
    let mut start = 0usize;
    let mut cell = 0usize;
    while start < nnz as usize {
        let count = 1_000_000usize.min(nnz as usize - start);
        let indices = read_integer_range(matrix.i, source, start, count)?;
        let values = read_real_range(matrix.x, source, start, count)?;
        for (offset, (&gene, &value)) in indices.iter().zip(&values).enumerate() {
            let position = start + offset;
            while cell + 1 < p.len() && position >= p[cell + 1] as usize {
                cell += 1;
            }
            ensure!(
                gene >= 0 && (gene as usize) < matrix.rows,
                "dgCMatrix row index {gene} is out of bounds"
            );
            if let Some(group) = cell_groups[cell] {
                sums[gene as usize * group_count + group] += to_linear(value, scale, matrix_name);
            }
        }
        start += count;
    }
    finish_group_means(&mut sums, &group_sizes);
    Ok(sums)
}

fn to_linear(value: f64, scale: InputScale, name: &str) -> f64 {
    match scale {
        InputScale::Linear => value,
        InputScale::Log1p => value.exp_m1(),
        InputScale::Auto if name.eq_ignore_ascii_case("counts") => value,
        InputScale::Auto => value.exp_m1(),
    }
}

fn seurat_metadata_frame(object: &RObject) -> Result<&rds2rust::DataFrameData> {
    ensure!(is_seurat(object), "RDS object is not a Seurat object");
    match seurat_slot(object, "meta.data")? {
        RObject::DataFrame(frame) => Ok(frame),
        other => bail!(
            "Seurat meta.data decoded as {}, not a data.frame",
            describe(other)
        ),
    }
}

fn seurat_slot<'a>(object: &'a RObject, name: &str) -> Result<&'a RObject> {
    as_s4(object, "root Seurat object")?
        .slots
        .get(name)
        .ok_or_else(|| anyhow!("Seurat object has no '{name}' slot"))
}

fn as_s4<'a>(object: &'a RObject, context: &str) -> Result<&'a rds2rust::S4ObjectData> {
    match object {
        RObject::S4Object(s4) => Ok(s4),
        RObject::WithAttributes { object, .. } => as_s4(object, context),
        other => bail!("{context} decoded as {}, not an S4 object", describe(other)),
    }
}

fn active_assay_name(object: &RObject) -> Result<String> {
    let value = seurat_slot(object, "active.assay")?;
    let values = loaded_characters(value, "active.assay")?;
    values
        .first()
        .and_then(|x| x.as_deref())
        .map(str::to_owned)
        .ok_or_else(|| anyhow!("active.assay is empty"))
}

fn seurat_assay<'a>(object: &'a RObject, name: &str) -> Result<&'a rds2rust::S4ObjectData> {
    let assay = named_item(seurat_slot(object, "assays")?, name)
        .with_context(|| format!("assay '{name}' does not exist"))?;
    as_s4(assay, &format!("assay '{name}'"))
}

fn seurat_assay_layout(assay: &rds2rust::S4ObjectData) -> Result<SeuratAssayLayout> {
    if assay
        .class
        .iter()
        .any(|class| matches!(class.as_ref(), "Assay5" | "StdAssay"))
    {
        Ok(SeuratAssayLayout::V5)
    } else if assay.class.iter().any(|class| class.as_ref() == "Assay") {
        Ok(SeuratAssayLayout::Legacy)
    } else {
        bail!(
            "assay class '{}' is unsupported; expected Assay or Assay5",
            join_classes(&assay.class)
        )
    }
}

fn seurat_format_summary(object: &RObject) -> Result<String> {
    ensure!(is_seurat(object), "RDS object is not a Seurat object");
    let active_assay = active_assay_name(object)?;
    let assay = seurat_assay(object, &active_assay)?;
    let layout = seurat_assay_layout(assay)?;
    let version = seurat_slot(object, "version")
        .ok()
        .and_then(render_package_version)
        .unwrap_or_else(|| "unknown version".to_owned());
    Ok(format!(
        "detected: Seurat {version}; active assay '{active_assay}' uses {}",
        layout.description()
    ))
}

fn render_package_version(object: &RObject) -> Option<String> {
    match object {
        RObject::Integer(VectorData::Owned(values)) if !values.is_empty() => {
            values.iter().all(|value| *value >= 0).then(|| {
                values
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(".")
            })
        }
        RObject::Character(VectorData::Owned(values)) => values
            .first()
            .and_then(|value| value.as_deref())
            .map(str::to_owned),
        RObject::List(values) if values.len() == 1 => render_package_version(&values[0]),
        RObject::S3Object(s3) => render_package_version(&s3.base),
        RObject::WithAttributes { object, .. } => render_package_version(object),
        _ => None,
    }
}

fn named_item<'a>(object: &'a RObject, wanted: &str) -> Result<&'a RObject> {
    let (value, attrs) = with_attributes(object)?;
    let items = match value {
        RObject::List(items) => items,
        other => bail!(
            "named collection decoded as {}, not a list",
            describe(other)
        ),
    };
    let names = loaded_characters(
        attrs
            .get("names")
            .ok_or_else(|| anyhow!("named collection has no names attribute"))?,
        "names",
    )?;
    let index = names
        .iter()
        .position(|name| name.as_deref() == Some(wanted))
        .ok_or_else(|| anyhow!("'{wanted}' is not present"))?;
    items
        .get(index)
        .ok_or_else(|| anyhow!("name index {index} is missing from collection"))
}

fn with_attributes(object: &RObject) -> Result<(&RObject, &Attributes)> {
    match object {
        RObject::WithAttributes { object, attributes } => Ok((object, attributes)),
        other => bail!("expected an attributed object, found {}", describe(other)),
    }
}

fn loaded_characters<'a>(object: &'a RObject, context: &str) -> Result<&'a [Option<Arc<str>>]> {
    match object {
        RObject::Character(VectorData::Owned(values)) => Ok(values),
        RObject::WithAttributes { object, .. } => loaded_characters(object, context),
        RObject::Character(VectorData::Lazy(_)) => bail!("{context} was unexpectedly left lazy"),
        other => bail!("{context} decoded as {}, not character", describe(other)),
    }
}

fn loaded_integers<'a>(object: &'a RObject, context: &str) -> Result<&'a [i32]> {
    match object {
        RObject::Integer(VectorData::Owned(values)) => Ok(values),
        RObject::WithAttributes { object, .. } => loaded_integers(object, context),
        RObject::Integer(VectorData::Lazy(_)) => bail!("{context} was unexpectedly left lazy"),
        other => bail!("{context} decoded as {}, not integer", describe(other)),
    }
}

fn sparse_matrix(object: &RObject) -> Result<SparseMatrixRef<'_>> {
    let s4 = as_s4(object, "expression layer")?;
    ensure!(
        s4.class.iter().any(|class| class.as_ref() == "dgCMatrix"),
        "layer class '{}' is unsupported; expected dgCMatrix",
        join_classes(&s4.class)
    );
    let dims = loaded_integers(
        s4.slots
            .get("Dim")
            .ok_or_else(|| anyhow!("dgCMatrix has no Dim slot"))?,
        "dgCMatrix Dim",
    )?;
    ensure!(
        dims.len() == 2 && dims[0] >= 0 && dims[1] >= 0,
        "invalid dgCMatrix dimensions"
    );
    let i = match s4.slots.get("i") {
        Some(RObject::Integer(values)) => values,
        _ => bail!("dgCMatrix has no integer i slot"),
    };
    let p = match s4.slots.get("p") {
        Some(RObject::Integer(values)) => values,
        _ => bail!("dgCMatrix has no integer p slot"),
    };
    let x = match s4.slots.get("x") {
        Some(RObject::Real(values)) => values,
        _ => bail!("dgCMatrix has no double x slot"),
    };
    Ok(SparseMatrixRef {
        rows: dims[0] as usize,
        cols: dims[1] as usize,
        i,
        p,
        x,
    })
}

fn layer_member_names(
    assay: &rds2rust::S4ObjectData,
    map_slot: &str,
    layer: &str,
    source: &dyn RdsInput,
) -> Result<Vec<String>> {
    let map = assay
        .slots
        .get(map_slot)
        .ok_or_else(|| anyhow!("Assay5 has no {map_slot} map"))?;
    let (values, attrs) = with_attributes(map)?;
    let logicals = match values {
        RObject::Logical(values) => values,
        other => bail!(
            "Assay5 {map_slot} map decoded as {}, not logical",
            describe(other)
        ),
    };
    let dims = loaded_integers(
        attrs
            .get("dim")
            .ok_or_else(|| anyhow!("{map_slot} map has no dimensions"))?,
        &format!("{map_slot} map dimensions"),
    )?;
    ensure!(
        dims.len() == 2 && dims[0] >= 0 && dims[1] >= 0,
        "invalid {map_slot} map dimensions"
    );
    let dimnames = match attrs.get("dimnames") {
        Some(RObject::List(values)) => values,
        _ => bail!("{map_slot} map has no dimnames"),
    };
    ensure!(dimnames.len() == 2, "invalid {map_slot} map dimnames");
    let all_names = read_all_characters(&dimnames[0], source)?;
    let layer_names = loaded_characters(&dimnames[1], &format!("{map_slot} layer names"))?;
    let layer_index = layer_names
        .iter()
        .position(|name| name.as_deref() == Some(layer))
        .ok_or_else(|| anyhow!("layer '{layer}' has no {map_slot} membership map"))?;
    let row_count = dims[0] as usize;
    ensure!(
        all_names.len() == row_count,
        "{map_slot} name count does not match map rows"
    );
    let membership = read_logical_range(logicals, source, layer_index * row_count, row_count)?;
    Ok(all_names
        .into_iter()
        .zip(membership)
        .filter_map(|(name, present)| match present {
            Logical::True => Some(name.as_deref().unwrap_or("NA").to_owned()),
            _ => None,
        })
        .collect())
}

fn seurat_cell_names(
    object: &RObject,
    source: &dyn RdsInput,
    count: usize,
) -> Result<Vec<Option<Arc<str>>>> {
    let active = active_assay_name(object)?;
    let assay = seurat_assay(object, &active)?;
    match seurat_assay_layout(assay)? {
        SeuratAssayLayout::V5 => {
            let cells = assay
                .slots
                .get("cells")
                .ok_or_else(|| anyhow!("Assay5 has no cells map"))?;
            let (_, attrs) = with_attributes(cells)?;
            let dimnames = match attrs.get("dimnames") {
                Some(RObject::List(values)) => values,
                _ => bail!("Assay5 cells map has no dimnames"),
            };
            read_character_range(&dimnames[0], source, 0, count)
        }
        SeuratAssayLayout::Legacy => {
            let matrix = assay
                .slots
                .get("data")
                .ok_or_else(|| anyhow!("legacy Assay has no data slot"))?;
            Ok(r_matrix_names(matrix, source, 1)?
                .into_iter()
                .take(count)
                .map(|name| Some(Arc::from(name)))
                .collect())
        }
    }
}

fn read_groups(
    object: &RObject,
    source: &dyn RdsInput,
) -> Result<(Vec<String>, Vec<Option<usize>>)> {
    if let Some((codes, levels)) = read_factor(object, source, None)? {
        let names = levels
            .iter()
            .map(|value| value.as_deref().unwrap_or("NA").to_owned())
            .collect::<Vec<_>>();
        let groups = codes
            .into_iter()
            .map(|code| {
                if code > 0 && code as usize <= names.len() {
                    Some(code as usize - 1)
                } else {
                    None
                }
            })
            .collect();
        return Ok((names, groups));
    }
    if matches!(object, RObject::Character(_)) {
        let values = read_all_characters(object, source)?;
        let mut names = Vec::<String>::new();
        let mut indices = HashMap::<String, usize>::new();
        let groups = values
            .into_iter()
            .map(|value| {
                value.map(|value| {
                    let key = value.to_string();
                    if let Some(index) = indices.get(&key) {
                        *index
                    } else {
                        let index = names.len();
                        names.push(key.clone());
                        indices.insert(key, index);
                        index
                    }
                })
            })
            .collect();
        return Ok((names, groups));
    }
    bail!(
        "grouping column must be a factor or character vector, found {}",
        describe(object)
    )
}

type FactorRead<'a> = Option<(Vec<i32>, &'a [Option<Arc<str>>])>;

fn read_factor<'a>(
    object: &'a RObject,
    source: &dyn RdsInput,
    count: Option<usize>,
) -> Result<FactorRead<'a>> {
    match object {
        RObject::Factor(factor) => {
            let take = count
                .unwrap_or(factor.values.len())
                .min(factor.values.len());
            Ok(Some((factor.values[..take].to_vec(), &factor.levels)))
        }
        RObject::S3Object(s3) if s3.class.iter().any(|class| class.as_ref() == "factor") => {
            let values = match s3.base.as_ref() {
                RObject::Integer(values) => values,
                _ => bail!("factor base is not integer"),
            };
            let levels = loaded_characters(
                s3.attributes
                    .get("levels")
                    .ok_or_else(|| anyhow!("factor has no levels"))?,
                "factor levels",
            )?;
            let take = count.unwrap_or(values.len()).min(values.len());
            Ok(Some((read_integer_range(values, source, 0, take)?, levels)))
        }
        _ => Ok(None),
    }
}

fn render_column_range(
    object: &RObject,
    source: &dyn RdsInput,
    count: usize,
) -> Result<Vec<String>> {
    if let Some((values, levels)) = read_factor(object, source, Some(count))? {
        return Ok(values
            .into_iter()
            .map(|code| {
                if code > 0 {
                    levels
                        .get(code as usize - 1)
                        .and_then(|x| x.as_deref())
                        .unwrap_or("NA")
                        .to_owned()
                } else {
                    "NA".into()
                }
            })
            .collect());
    }
    match object {
        RObject::Character(_) => Ok(read_character_range(object, source, 0, count)?
            .into_iter()
            .map(|x| x.as_deref().unwrap_or("NA").to_owned())
            .collect()),
        RObject::Integer(values) => Ok(read_integer_range(values, source, 0, count)?
            .into_iter()
            .map(|x| {
                if x == i32::MIN {
                    "NA".into()
                } else {
                    x.to_string()
                }
            })
            .collect()),
        RObject::Real(values) => Ok(read_real_range(values, source, 0, count)?
            .into_iter()
            .map(|x| {
                if x.is_nan() {
                    "NA".into()
                } else {
                    x.to_string()
                }
            })
            .collect()),
        RObject::Logical(values) => Ok(read_logical_range(values, source, 0, count)?
            .into_iter()
            .map(|x| {
                match x {
                    Logical::True => "TRUE",
                    Logical::False => "FALSE",
                    Logical::Na => "NA",
                }
                .into()
            })
            .collect()),
        other => Ok(vec![format!("<{}>", describe(other)); count]),
    }
}

fn read_integer_range(
    data: &VectorData<i32>,
    source: &dyn RdsInput,
    start: usize,
    count: usize,
) -> Result<Vec<i32>> {
    match data {
        VectorData::Owned(values) => Ok(values[start..start + count].to_vec()),
        VectorData::Lazy(span) => Ok(read_lazy_integer_range(source, *span, start, count)?),
    }
}

fn read_real_range(
    data: &VectorData<f64>,
    source: &dyn RdsInput,
    start: usize,
    count: usize,
) -> Result<Vec<f64>> {
    match data {
        VectorData::Owned(values) => Ok(values[start..start + count].to_vec()),
        VectorData::Lazy(span) => Ok(read_lazy_real_range(source, *span, start, count)?),
    }
}

fn read_logical_range(
    data: &VectorData<Logical>,
    source: &dyn RdsInput,
    start: usize,
    count: usize,
) -> Result<Vec<Logical>> {
    match data {
        VectorData::Owned(values) => Ok(values[start..start + count].to_vec()),
        VectorData::Lazy(span) => Ok(read_lazy_logical_range(source, *span, start, count)?),
    }
}

fn read_character_range(
    object: &RObject,
    source: &dyn RdsInput,
    start: usize,
    count: usize,
) -> Result<Vec<Option<Arc<str>>>> {
    let data = match object {
        RObject::Character(data) => data,
        RObject::WithAttributes { object, .. } => {
            return read_character_range(object, source, start, count);
        }
        other => bail!("expected character vector, found {}", describe(other)),
    };
    match data {
        VectorData::Owned(values) => Ok(values[start..start + count].to_vec()),
        VectorData::Lazy(span) => Ok(read_lazy_character_range(source, *span, start, count)?),
    }
}

fn read_all_integers(data: &VectorData<i32>, source: &dyn RdsInput) -> Result<Vec<i32>> {
    read_integer_range(data, source, 0, data.len())
}

fn read_all_characters(object: &RObject, source: &dyn RdsInput) -> Result<Vec<Option<Arc<str>>>> {
    let length = match object {
        RObject::Character(data) => data.len(),
        RObject::WithAttributes { object, .. } => return read_all_characters(object, source),
        other => bail!("expected character vector, found {}", describe(other)),
    };
    read_character_range(object, source, 0, length)
}

fn default_output(input: &Path) -> PathBuf {
    let mut output = input.with_extension("");
    let name = output
        .file_name()
        .and_then(|x| x.to_str())
        .unwrap_or("output");
    output.set_file_name(format!("{name}.refmat.tsv"));
    output
}

fn write_reference(
    path: &Path,
    features: &[String],
    groups: &[String],
    values: &[f64],
) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    write!(writer, "gene")?;
    for group in groups {
        write!(writer, "\t{}", escape_tsv(group))?;
    }
    writeln!(writer)?;
    for (gene, feature) in features.iter().enumerate() {
        write!(writer, "{}", escape_tsv(feature))?;
        for group in 0..groups.len() {
            write!(writer, "\t{}", values[gene * groups.len() + group])?;
        }
        writeln!(writer)?;
    }
    writer.flush()?;
    Ok(())
}

fn escape_tsv(value: &str) -> String {
    value.replace(['\t', '\n', '\r'], " ")
}

fn print_table<'a>(
    row_header: &str,
    row_names: &[String],
    column_names: impl IntoIterator<Item = &'a str>,
    columns: &[Vec<String>],
) {
    let column_names = column_names.into_iter().collect::<Vec<_>>();
    debug_assert_eq!(column_names.len(), columns.len());
    debug_assert!(columns.iter().all(|column| column.len() == row_names.len()));

    let mut widths = Vec::with_capacity(columns.len() + 1);
    widths.push(
        row_names
            .iter()
            .map(|value| display_width(value))
            .max()
            .unwrap_or(0)
            .max(display_width(row_header)),
    );
    widths.extend(column_names.iter().zip(columns).map(|(name, values)| {
        values
            .iter()
            .map(|value| display_width(value))
            .max()
            .unwrap_or(0)
            .max(display_width(name))
    }));

    let mut output = String::new();
    write_table_rule(&mut output, &widths);
    write_table_row(
        &mut output,
        std::iter::once(row_header).chain(column_names.iter().copied()),
        &widths,
    );
    write_table_rule(&mut output, &widths);
    for row in 0..row_names.len() {
        write_table_row(
            &mut output,
            std::iter::once(row_names[row].as_str())
                .chain(columns.iter().map(|column| column[row].as_str())),
            &widths,
        );
    }
    write_table_rule(&mut output, &widths);
    print!("{output}");
}

const METADATA_TABLE_WIDTH: usize = 120;

fn print_metadata_table<'a>(
    row_names: &[String],
    column_names: impl IntoIterator<Item = &'a str>,
    columns: &[Vec<String>],
) {
    let column_names = column_names.into_iter().collect::<Vec<_>>();
    if columns.is_empty() {
        print_table("cell", row_names, std::iter::empty(), columns);
        return;
    }
    let ranges = metadata_column_ranges(row_names, &column_names, columns, METADATA_TABLE_WIDTH);
    let split = ranges.len() > 1;

    for (block, range) in ranges.into_iter().enumerate() {
        if block > 0 {
            println!();
        }
        if split {
            println!(
                "Metadata columns {}-{} of {}:",
                range.start + 1,
                range.end,
                columns.len()
            );
        }
        print_table(
            "cell",
            row_names,
            column_names[range.clone()].iter().copied(),
            &columns[range],
        );
    }
}

fn metadata_column_ranges(
    row_names: &[String],
    column_names: &[&str],
    columns: &[Vec<String>],
    max_width: usize,
) -> Vec<std::ops::Range<usize>> {
    debug_assert_eq!(column_names.len(), columns.len());
    if columns.is_empty() {
        return Vec::new();
    }

    let row_width = row_names
        .iter()
        .map(|value| display_width(value))
        .max()
        .unwrap_or(0)
        .max("cell".len());
    let mut ranges = Vec::new();
    let mut start = 0;
    while start < columns.len() {
        let mut end = start;
        let mut width = row_width + 4;
        while end < columns.len() {
            let column_width = columns[end]
                .iter()
                .map(|value| display_width(value))
                .max()
                .unwrap_or(0)
                .max(display_width(column_names[end]));
            let added = column_width + 3;
            if end > start && width + added > max_width {
                break;
            }
            width += added;
            end += 1;
        }
        ranges.push(start..end);
        start = end;
    }
    ranges
}

fn print_group_counts(column: &str, group_names: &[String], cell_groups: &[Option<usize>]) {
    let mut counts = vec![0usize; group_names.len()];
    for group in cell_groups.iter().flatten() {
        if let Some(count) = counts.get_mut(*group) {
            *count += 1;
        }
    }
    let counts = counts
        .into_iter()
        .map(|count| count.to_string())
        .collect::<Vec<_>>();
    print_table(column, group_names, ["cells"], &[counts]);
}

fn display_width(value: &str) -> usize {
    escape_tsv(value).chars().count()
}

fn write_table_rule(output: &mut String, widths: &[usize]) {
    output.push('+');
    for width in widths {
        output.push_str(&"-".repeat(width + 2));
        output.push('+');
    }
    output.push('\n');
}

fn write_table_row<'a>(
    output: &mut String,
    values: impl IntoIterator<Item = &'a str>,
    widths: &[usize],
) {
    output.push('|');
    for (value, width) in values.into_iter().zip(widths) {
        let value = escape_tsv(value);
        let padding = width.saturating_sub(value.chars().count());
        let _ = write!(output, " {value}{} |", " ".repeat(padding));
    }
    output.push('\n');
}

fn finish_group_means(sums: &mut [f64], group_sizes: &[usize]) {
    for row in sums.chunks_exact_mut(group_sizes.len()) {
        for (value, size) in row.iter_mut().zip(group_sizes) {
            *value = (*value / *size as f64).ln_1p();
        }
    }
}

fn print_object(label: &str, object: &RObject, level: usize, max_depth: usize) {
    let indent = "  ".repeat(level);
    println!("{indent}{label}: {}", describe(object));
    if level >= max_depth {
        return;
    }

    match object {
        RObject::S4Object(s4) => {
            for (name, value) in &s4.slots {
                print_object(&format!("@{name}"), value, level + 1, max_depth);
            }
        }
        RObject::S3Object(s3) => {
            print_object("base", &s3.base, level + 1, max_depth);
            for (name, value) in s3.attributes.iter() {
                print_object(&format!("attr({name})"), value, level + 1, max_depth);
            }
        }
        RObject::DataFrame(frame) => {
            for (name, value) in &frame.columns {
                print_object(name, value, level + 1, max_depth);
            }
        }
        RObject::List(values) | RObject::Expression(values) => {
            for (index, value) in values.iter().enumerate().take(50) {
                print_object(&format!("[[{}]]", index + 1), value, level + 1, max_depth);
            }
            if values.len() > 50 {
                println!("{}  ... {} more elements", indent, values.len() - 50);
            }
        }
        RObject::WithAttributes { object, attributes } => {
            print_object("value", object, level + 1, max_depth);
            for (name, value) in attributes.iter() {
                print_object(&format!("attr({name})"), value, level + 1, max_depth);
            }
        }
        RObject::Shared(shared) => match shared.read() {
            Ok(value) => print_object("shared", &value, level + 1, max_depth),
            Err(_) => println!("{}  <poisoned shared reference>", indent),
        },
        _ => {}
    }
}

fn describe(object: &RObject) -> String {
    match object {
        RObject::Null => "NULL".into(),
        RObject::Integer(v) => describe_vector("integer", v),
        RObject::Real(v) => describe_vector("double", v),
        RObject::Logical(v) => describe_vector("logical", v),
        RObject::Character(v) => describe_vector("character", v),
        RObject::Raw(v) => describe_vector("raw", v),
        RObject::Complex(v) => describe_vector("complex", v),
        RObject::Symbol(name) => format!("symbol({name})"),
        RObject::List(values) => format!("list[{}]", values.len()),
        RObject::Pairlist(values) => format!("pairlist[{}]", values.len()),
        RObject::Expression(values) => format!("expression[{}]", values.len()),
        RObject::DataFrame(frame) => format!(
            "data.frame[{} rows x {} columns]",
            frame.row_names.len(),
            frame.columns.len()
        ),
        RObject::Factor(factor) => format!(
            "{}factor[{}; {} levels]",
            if factor.ordered { "ordered " } else { "" },
            factor.values.len(),
            factor.levels.len()
        ),
        RObject::S3Object(s3) => format!("S3<{}>", join_classes(&s3.class)),
        RObject::S4Object(s4) => format!(
            "S4<{}>{}",
            join_classes(&s4.class),
            s4.package
                .as_deref()
                .map_or_else(String::new, |p| format!(" from {p}"))
        ),
        RObject::Shared(_) => "shared reference".into(),
        RObject::WithAttributes { attributes, .. } => {
            let names = attributes
                .iter()
                .map(|(name, _)| name.as_ref())
                .collect::<Vec<_>>()
                .join(", ");
            format!("object with attributes [{names}]")
        }
        RObject::Environment { .. } => "environment".into(),
        RObject::Closure { .. } => "closure".into(),
        RObject::Language { .. } => "language".into(),
        RObject::Promise { .. } => "promise".into(),
        RObject::Special { name } => format!("special({name})"),
        RObject::Builtin { name } => format!("builtin({name})"),
        RObject::Bytecode { .. } => "bytecode".into(),
        RObject::Namespace(name) => format!("namespace({})", join_classes(name)),
        RObject::PackageEnv(name) => format!("package-env({})", join_classes(name)),
        RObject::GlobalEnv => "global environment".into(),
        RObject::BaseEnv => "base environment".into(),
        RObject::EmptyEnv => "empty environment".into(),
        RObject::MissingArg => "missing argument".into(),
        RObject::UnboundValue => "unbound value".into(),
        _ => "unsupported/unknown R object".into(),
    }
}

fn describe_vector<T>(kind: &str, data: &VectorData<T>) -> String {
    let storage = if data.is_loaded() { "loaded" } else { "lazy" };
    format!("{kind}[{}] ({storage})", data.len())
}

fn join_classes(values: &[Arc<str>]) -> String {
    values
        .iter()
        .map(AsRef::as_ref)
        .collect::<Vec<_>>()
        .join(", ")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finishes_linear_group_sums_on_log_scale() {
        let mut sums = vec![3.0, 8.0, 0.0, 2.0];
        finish_group_means(&mut sums, &[3, 2]);
        let expected = [1.0_f64.ln_1p(), 4.0_f64.ln_1p(), 0.0, 1.0_f64.ln_1p()];
        for (actual, expected) in sums.iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-15);
        }
    }

    #[test]
    fn derives_default_output_next_to_input() {
        assert_eq!(
            default_output(Path::new("data/so.rds")),
            PathBuf::from("data/so.refmat.tsv")
        );
    }

    #[test]
    fn escapes_tsv_control_characters() {
        assert_eq!(escape_tsv("a\tb\nc\r"), "a b c ");
    }

    #[test]
    fn splits_wide_metadata_into_readable_column_blocks() {
        let rows = vec!["Cell1".to_owned()];
        let names = ["one", "two", "three"];
        let columns = vec![
            vec!["1234".to_owned()],
            vec!["1234".to_owned()],
            vec!["1234".to_owned()],
        ];
        assert_eq!(
            metadata_column_ranges(&rows, &names, &columns, 25),
            vec![0..2, 2..3]
        );
    }
}
