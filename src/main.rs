use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, anyhow, bail, ensure};
use clap::{Parser, Subcommand};
use rds2rust::{
    Attributes, ChunkedRdsSource, Logical, ParseConfig, RObject, RdsInput, VectorData,
    read_lazy_character_range, read_lazy_integer_range, read_lazy_logical_range,
    read_lazy_real_range, read_rds_from_path_chunked, read_rds_with_input,
};

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
    /// Build a genes-by-group reference matrix.
    Build {
        file: PathBuf,
        #[arg(short, long)]
        column: String,
        #[arg(long)]
        assay: Option<String>,
        #[arg(long)]
        layer: Option<String>,
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Inspect { file, depth, full } => inspect(&file, depth, full),
        Command::Head { file, rows } => head(&file, rows),
        Command::Build {
            file,
            column,
            assay,
            layer,
            output,
        } => build(
            &file,
            &column,
            assay.as_deref(),
            layer.as_deref(),
            output.as_deref(),
        ),
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
    if full {
        let object = read_rds_from_path_chunked(file)
            .with_context(|| format!("failed to parse {}", file.display()))?
            .object;
        print_object("root", &object, 0, depth);
    } else {
        let parsed = parse_lazy(file)?;
        eprintln!("note: {} large vectors kept lazy", parsed.lazy_vectors);
        print_object("root", &parsed.object, 0, depth);
    }
    Ok(())
}

fn head(file: &Path, rows: usize) -> Result<()> {
    let parsed = parse_lazy(file)?;
    let frame = metadata_frame(&parsed.object)?;
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

    print!("cell");
    for name in frame.columns.keys() {
        print!("\t{name}");
    }
    println!();
    for row in 0..take {
        print!("{}", escape_tsv(row_names[row].as_deref().unwrap_or("NA")));
        for column in &rendered {
            print!("\t{}", escape_tsv(&column[row]));
        }
        println!();
    }
    Ok(())
}

fn build(
    file: &Path,
    column: &str,
    assay: Option<&str>,
    layer: Option<&str>,
    output: Option<&Path>,
) -> Result<()> {
    let parsed = parse_lazy(file)?;
    let frame = metadata_frame(&parsed.object)?;
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
    let assay_object = named_item(seurat_slot(&parsed.object, "assays")?, assay_name)
        .with_context(|| format!("assay '{assay_name}' does not exist"))?;
    let assay_s4 = as_s4(assay_object, "assay")?;
    let layer_object = named_item(
        assay_s4
            .slots
            .get("layers")
            .ok_or_else(|| anyhow!("Assay5 has no layers slot"))?,
        layer_name,
    )
    .with_context(|| format!("layer '{layer_name}' does not exist in assay '{assay_name}'"))?;
    let matrix = sparse_matrix(layer_object)?;

    let cell_names = layer_member_names(assay_s4, "cells", layer_name, &parsed.source)?;
    let feature_names = layer_member_names(assay_s4, "features", layer_name, &parsed.source)?;
    ensure!(
        matrix.cols == cell_names.len(),
        "layer has {} columns but its cell map selects {} cells",
        matrix.cols,
        cell_names.len()
    );
    ensure!(
        matrix.rows == feature_names.len(),
        "layer has {} rows but its feature map selects {} features",
        matrix.rows,
        feature_names.len()
    );
    ensure!(
        cell_groups.len() == matrix.cols,
        "metadata has {} cells but layer has {}; partial-layer cell alignment is not implemented yet",
        cell_groups.len(),
        matrix.cols
    );

    let mut group_sizes = vec![0usize; group_names.len()];
    for group in cell_groups.iter().flatten() {
        group_sizes[*group] += 1;
    }
    ensure!(
        group_sizes.iter().all(|size| *size > 0),
        "annotation contains an empty factor level"
    );

    let p = read_all_integers(matrix.p, &parsed.source)?;
    ensure!(
        p.len() == matrix.cols + 1,
        "invalid dgCMatrix p length {} for {} columns",
        p.len(),
        matrix.cols
    );
    ensure!(p.first() == Some(&0), "invalid dgCMatrix: p[0] is not zero");
    ensure!(
        p.windows(2).all(|x| x[0] <= x[1]),
        "invalid dgCMatrix: p is not monotonic"
    );
    let nnz = *p.last().unwrap_or(&0);
    ensure!(
        nnz >= 0 && nnz as usize == matrix.i.len() && nnz as usize == matrix.x.len(),
        "invalid dgCMatrix nonzero lengths"
    );

    let mut sums = vec![0.0f64; matrix.rows * group_names.len()];
    let chunk_size = 1_000_000usize;
    let mut start = 0usize;
    let mut cell = 0usize;
    while start < nnz as usize {
        let count = chunk_size.min(nnz as usize - start);
        let indices = read_integer_range(matrix.i, &parsed.source, start, count)?;
        let values = read_real_range(matrix.x, &parsed.source, start, count)?;
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
                let linear = if layer_name == "counts" {
                    value
                } else {
                    value.exp_m1()
                };
                sums[gene as usize * group_names.len() + group] += linear;
            }
        }
        start += count;
    }

    finish_group_means(&mut sums, &group_sizes);

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

fn metadata_frame(object: &RObject) -> Result<&rds2rust::DataFrameData> {
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
    let assay = named_item(seurat_slot(object, "assays")?, &active)?;
    let assay = as_s4(assay, "active assay")?;
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
}
